% FILE: run_voxel_arima_jhpce.m
function run_voxel_arima_jhpce()
    % This script is designed to be run non-interactively on the JHPCE cluster.
    % It saves all results to a .mat file and NIfTI files.

    fprintf('Starting ARIMA fitting job at %s\n', datestr(now));

    %% --- 0. Cluster-Specific Setup for Parallel Pool ---
    % Automatically detect the number of cores allocated by SLURM
    c = parcluster;
    num_cores_str = getenv('SLURM_CPUS_PER_TASK');
    if isempty(num_cores_str)
        warning('SLURM_CPUS_PER_TASK not set. Using default parallel pool settings.');
        num_cores = c.NumWorkers; % Fallback for local testing
    else
        num_cores = str2double(num_cores_str);
    end
    fprintf('Requesting a parallel pool with %d workers.\n', num_cores);
    % Start the parallel pool with the correct number of workers
    parpool(c, num_cores);


    %% --- 1. Load Data and Identify Target Voxels ---
    % Add CanlabCore to path (adjust if your startup.m doesn't do this)
    % addpath(genpath('/path/to/your/CanlabCore'));

    % Load the 4D residual data file
    dat = fmri_data('data/WM_100307_res4d.nii');
    time_series_data = double(dat.dat);

    % --- This part is time-consuming, so let's check if we've done it before ---
    lbq_results_file = 'lbq_results.mat';
    if exist(lbq_results_file, 'file')
        fprintf('Loading pre-computed Ljung-Box results...\n');
        load(lbq_results_file, 'lbq_map');
    else
        fprintf('Running Ljung-Box tests to find significant voxels...\n');
        % (Insert the Ljung-Box test code from our previous conversation here)
        % ...
        % At the end of that code, you should have the variable 'lbq_map'
        fprintf('Saving Ljung-Box results to %s\n', lbq_results_file);
        save(lbq_results_file, 'lbq_map');
    end

    alpha = 0.05;
    n_voxels_total = size(time_series_data, 1);
    sig_voxels_idx = find(lbq_map < alpha);
    n_sig_voxels = length(sig_voxels_idx);

    if n_sig_voxels == 0
        fprintf('No voxels found with significant autocorrelation. Exiting.\n');
        return;
    end
    fprintf('Found %d voxels with significant autocorrelation to model.\n', n_sig_voxels);

    valid_residuals = detrend(time_series_data(sig_voxels_idx, :)')';

    %% --- 2. Setup for ARIMA Fitting ---
    s = 5;      max_p = 2;  max_q = 2;  max_P = 1;  max_Q = 1;

    %% --- 3. Main Loop: Fit Models to Significant Voxels ---
    best_p_vals = nan(n_sig_voxels, 1);
    best_q_vals = nan(n_sig_voxels, 1);
    best_P_vals = nan(n_sig_voxels, 1);
    best_Q_vals = nan(n_sig_voxels, 1);
    best_bic_vals = nan(n_sig_voxels, 1);

    fprintf('Starting main parallel loop for ARIMA fitting...\n');
    tic;
    parfor i = 1:n_sig_voxels
        % The parfor loop itself is quiet, which is fine for batch jobs.
        voxel_ts = valid_residuals(i, :);
        [~, p, q, P, Q, bic] = find_best_sarima_model(voxel_ts, s, max_p, max_q, max_P, max_Q);

        best_p_vals(i) = p;
        best_q_vals(i) = q;
        best_P_vals(i) = P;
        best_Q_vals(i) = Q;
        best_bic_vals(i) = bic;
    end
    elapsed_time = toc;
    fprintf('ARIMA fitting completed in %.2f hours.\n', elapsed_time / 3600);

    %% --- 4. Reconstruct and Save Full Brain Images ---
    fprintf('Reconstructing and saving result maps...\n');
    p_map = nan(n_voxels_total, 1); q_map = nan(n_voxels_total, 1);
    P_map = nan(n_voxels_total, 1); Q_map = nan(n_voxels_total, 1);
    bic_map = nan(n_voxels_total, 1);

    p_map(sig_voxels_idx) = best_p_vals;   q_map(sig_voxels_idx) = best_q_vals;
    P_map(sig_voxels_idx) = best_P_vals;   Q_map(sig_voxels_idx) = best_Q_vals;
    bic_map(sig_voxels_idx) = best_bic_vals;
    
    % --- Save results to a .mat file for easy inspection ---
    save('arima_results.mat', 'p_map', 'q_map', 'P_map', 'Q_map', 'bic_map', 'sig_voxels_idx', 'elapsed_time');

    % --- Save results as NIfTI files for visualization ---
    p_obj = dat; p_obj.dat = p_map; p_obj.dat_descrip = 'Best AR order (p)';
    write(p_obj, 'fname', 'p_map.nii', 'overwrite');

    q_obj = dat; q_obj.dat = q_map; q_obj.dat_descrip = 'Best MA order (q)';
    write(q_obj, 'fname', 'q_map.nii', 'overwrite');
    
    P_obj = dat; P_obj.dat = P_map; P_obj.dat_descrip = 'Best SAR order (P)';
    write(P_obj, 'fname', 'P_map.nii', 'overwrite');

    Q_obj = dat; Q_obj.dat = Q_map; Q_obj.dat_descrip = 'Best SMA order (Q)';
    write(Q_obj, 'fname', 'Q_map.nii', 'overwrite');

    bic_obj = dat; bic_obj.dat = bic_map; bic_obj.dat_descrip = 'Best Model BIC';
    write(bic_obj, 'fname', 'bic_map.nii', 'overwrite');

    fprintf('Job finished successfully at %s\n', datestr(now));
    
    % Exit MATLAB
    exit;
end