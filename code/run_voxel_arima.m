%% --- 0. Initial Setup ---
addpath(genpath('/users/rwang/CanLabCore'));
addpath('/users/rwang/matlab_utilities'); 

fprintf('--- Starting Voxel-wise ARMA Analysis ---\n');

% --- Directory and subject setup ---
base_dir = '/dcs07/smart/data/HRF/bogdan_hrf/stats/fsl_noar_24vec_csf_spikes/results/WM';
output_dir = fullfile(base_dir, 'autocorrelation_analysis_output');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

all_items = dir(base_dir);
all_dirs = all_items([all_items.isdir]);
is_dot_dir = ismember({all_dirs.name}, {'.', '..'});
all_dirs = all_dirs(~is_dot_dir);

dir_names = {all_dirs.name};
is_numeric_6digit = ~cellfun('isempty', regexp(dir_names, '^\d{6}$'));
subject_dirs = all_dirs(is_numeric_6digit);

if isempty(subject_dirs)
    error('No 6-digit numeric subject directories found in %s.', base_dir);
end
subject_ids = {subject_dirs.name};
n_subjects = length(subject_ids);

% --- Get the subject index from the Slurm job array task ID ---
task_id_str = getenv('SLURM_ARRAY_TASK_ID');
if isempty(task_id_str)
    error('This script is designed for a Slurm job array. SLURM_ARRAY_TASK_ID not found.');
end
i = str2double(task_id_str);

if i > n_subjects
    fprintf('Task ID %d exceeds the number of subjects (%d). Exiting.\n', i, n_subjects);
    return;
end

subject_id = subject_ids{i};
fprintf('Processing Subject %d/%d: %s\n', i, n_subjects, subject_id);


%% --- 1. Load and Combine Data ---
file_LR = fullfile(base_dir, subject_id, 'LR', 'res4d.nii');
file_RL = fullfile(base_dir, subject_id, 'RL', 'res4d.nii');

if ~exist(file_LR, 'file') || ~exist(file_RL, 'file')
    fprintf('Data files for subject %s not found. Skipping.\n', subject_id);
    return;
end

fprintf('Loading data for subject %s...\n', subject_id);
data_LR = fmri_data(file_LR);
data_RL = fmri_data(file_RL);

combined_ts_data = [data_LR.dat, data_RL.dat];
[n_voxels, T] = size(combined_ts_data);
combined_ts_data_T = combined_ts_data';

fprintf('Data loaded. Found %d voxels and %d total time points.\n', n_voxels, T);


%% --- 2. Voxel-wise ARMA Model Fitting ---
max_p = 2;  
max_q = 2;  
p_vals = 0:max_p;
q_vals = 0:max_q;
ar_order_map = zeros(1, n_voxels);
ma_order_map = zeros(1, n_voxels);
error_messages = cell(n_voxels, 1);

if isempty(gcp('nocreate'))
    num_cores_str = getenv('SLURM_CPUS_PER_TASK');
    if isempty(num_cores_str)
        num_cores = feature('numcores'); % Fallback for local testing
        fprintf('Not on Slurm, starting parpool with %d workers.\n', num_cores);
    else
        num_cores = str2double(num_cores_str);
        fprintf('Detected %d CPUs from Slurm, starting parpool...\n', num_cores);
    end
    parpool(num_cores);
end
fprintf('Starting parallel ARMA model fitting for %d voxels...\n', n_voxels);

progress_file = fullfile(output_dir, [subject_id '_progress.txt']);
parfor_progress(n_voxels, progress_file);

parfor v = 1:n_voxels
    voxel_ts = double(combined_ts_data_T(:, v));
    voxel_ts = detrend(voxel_ts);
    if std(voxel_ts) > 1e-10 
        voxel_ts = voxel_ts / std(voxel_ts);
    end

    if var(voxel_ts) < 1e-4
        parfor_progress;
        continue;
    end

    best_bic = inf;
    best_p = 0;
    best_q = 0;

    for p = p_vals
        for q = q_vals
            try
                Mdl = arima(p, 0, q);
                [~, ~, logL] = estimate(Mdl, voxel_ts, 'Display', 'off', 'MaxIterations', 50);
                num_params = p + q + 1; % AR params + MA params + variance
                [~, bic] = aicbic(logL, num_params, T);

                if bic < best_bic
                    best_bic = bic;
                    best_p = p;
                    best_q = q;
                end

            catch ME
                error_messages{v} = sprintf('Voxel %d, p=%d, q=%d failed: %s', v, p, q, ME.message);
            end
        end
    end
    ar_order_map(v) = best_p;
    ma_order_map(v) = best_q;
    
    parfor_progress;
end

parfor_progress(0);

fprintf('Finished model fitting.\n');

% Log any captured errors
fprintf('\n--- ERROR LOG ---\n');
captured_errors = error_messages(~cellfun('isempty', error_messages));
num_errors = length(captured_errors);
if num_errors > 0
    max_errors_to_log = 10;
    fprintf('Captured %d errors. Displaying up to %d:\n', num_errors, max_errors_to_log);
    for err_idx = 1:min(num_errors, max_errors_to_log)
        fprintf('%s\n', captured_errors{err_idx});
    end
else
    fprintf('No errors were captured.\n');
end
fprintf('--- END ERROR LOG ---\n\n');


%% --- 3. Save Results as NIFTI Brain Maps ---
fprintf('Saving AR and MA order maps...\n');

template_obj = fmri_data(file_LR);

% Create and write the AR order map
ar_map_obj = statistic_image('dat', ar_order_map', 'volInfo', template_obj.volInfo, 'dat_descrip', ['AR order map for subject ' subject_id]);
ar_map_filename = fullfile(output_dir, [subject_id '_ar_order_map.nii']);
write(ar_map_obj, 'fname', ar_map_filename, 'overwrite');

% Create and write the MA order map
ma_map_obj = statistic_image('dat', ma_order_map', 'volInfo', template_obj.volInfo, 'dat_descrip', ['MA order map for subject ' subject_id]);
ma_map_filename = fullfile(output_dir, [subject_id '_ma_order_map.nii']);
write(ma_map_obj, 'fname', ma_map_filename, 'overwrite');

fprintf('--- Analysis for subject %s complete. Output maps saved in: %s ---\n', subject_id, output_dir);