%% --- 0. Initial Setup ---
% This script assumes that CanLabCore is on your MATLAB path.
% If not, add it now, e.g., addpath(genpath('C:\path\to\CanLabCore'));

% Define the directory containing your NIfTI files
data_dir = 'C:\Users\Ryan_\Documents\my_Git\fMRI\data';

% Create a directory for the output files to keep things organized
output_dir = fullfile(data_dir, 'autocorrelation_analysis_output');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Get a list of the residual files to be processed
subject_files = dir(fullfile(data_dir, 'WM_*_res4d.nii'));
n_subjects = length(subject_files);

if n_subjects == 0
    error('No residual files found. Please check the data_dir path and file names.');
end

% --- Start Parallel Pool ---
% Start a parallel pool for faster processing if one isn't already running.
if isempty(gcp('nocreate'))
    parpool;
end


%% --- 1. Subject-Level Analysis: Create Autocorrelation Maps ---

fprintf('--- Starting Step 1: Subject-Level Autocorrelation Mapping ---\n\n');

% Pre-allocate a cell array to store the paths to the generated Q-statistic maps
% This will be used in the group analysis step.
q_stat_filenames = cell(n_subjects, 1);

for i = 1:n_subjects
    subject_filepath = fullfile(subject_files(i).folder, subject_files(i).name);
    fprintf('Processing Subject %d/%d: %s\n', i, n_subjects, subject_files(i).name);

    % Load the 4D residual data for one subject using fmri_data
    dat = fmri_data(subject_filepath);
    time_series_data = double(dat.dat);

    % --- Setup for this subject's test ---
    [n_voxels_total, n_time_points] = size(time_series_data);
    % Set lags for the Ljung-Box test. A common heuristic is log(T) or up to 20.
    lags_for_lbq = min(20, floor(log(n_time_points)));
    fprintf('  Using %d lags for Ljung-Box test.\n', lags_for_lbq);

    % --- Data Cleaning ---
    % Find valid voxels (non-constant variance and no NaNs) to avoid errors.
    % This is a crucial step for robust processing.
    valid_voxels_mask = var(time_series_data, 0, 2) > 1e-6 & ~any(isnan(time_series_data), 2);
    data_valid = time_series_data(valid_voxels_mask, :);
    n_voxels_valid = size(data_valid, 1);
    
    % The residuals from a GLM should be mean-zero, but we detrend to be safe.
    preprocessed_residuals = detrend(data_valid', 'constant')';

    % --- Voxel-wise Ljung-Box Test (using a parallel for-loop) ---
    fprintf('  Running Ljung-Box test on %d valid voxels...\n', n_voxels_valid);
    
    % Pre-allocate result vectors for valid voxels
    q_stats_vec = nan(n_voxels_valid, 1);
    p_values_vec = nan(n_voxels_valid, 1);

    parfor v = 1:n_voxels_valid
        % For each voxel's time series...
        voxel_series = preprocessed_residuals(v, :);
        
        % Run the Ljung-Box Q-test. We are interested in the test statistic
        % and p-value for the maximum number of lags specified.
        [~, p, stat, ~] = lbqtest(voxel_series, 'Lags', lags_for_lbq);
        
        q_stats_vec(v) = stat(end);  % The Q-statistic for the max lag
        p_values_vec(v) = p(end);    % The p-value for the max lag
    end

    % --- Reconstruct Full Brain Maps ---
    % Create full-brain vectors, initializing with NaNs
    q_stat_map = nan(n_voxels_total, 1);
    p_value_map = nan(n_voxels_total, 1);
    
    % Place the calculated results back into their correct voxel locations
    q_stat_map(valid_voxels_mask) = q_stats_vec;
    p_value_map(valid_voxels_mask) = p_values_vec;

    % --- Create and Save fmri_data Objects ---
    [~, basename, ~] = fileparts(subject_files(i).name);
    
    % Create an fmri_data object for the Q-statistic map
    q_stat_obj = dat; % Use original object as a template for header info
    q_stat_obj.dat = q_stat_map;
    q_stat_obj.dat_descrip = 'Ljung-Box Q-statistic';
    
    % Define the output filename and save the NIfTI file
    output_q_filename = fullfile(output_dir, [basename '_Q_map.nii']);
    write(q_stat_obj, 'fname', output_q_filename, 'overwrite');
    q_stat_filenames{i} = output_q_filename; % Store filename for group analysis
    
    fprintf('  Saved Q-statistic map to: %s\n', output_q_filename);
    
    % (Optional) Save the p-value map as well for reference
    p_value_obj = dat;
    p_value_obj.dat = p_value_map;
    p_value_obj.dat_descrip = 'Ljung-Box p-value';
    output_p_filename = fullfile(output_dir, [basename '_p_map.nii']);
    write(p_value_obj, 'fname', output_p_filename, 'overwrite');
    fprintf('  Saved p-value map to: %s\n\n', output_p_filename);

end % End of subject loop

%% --- 2. Group-Level Analysis: One-Sample T-test ---

% --- Concatenate Individual Maps into a 4D object ---
fprintf('--- Starting Step 2: Group-Level Analysis ---\n\n');
fprintf('Concatenating individual Q-statistic maps...\n');
group_q_stats_data = fmri_data(q_stat_filenames);

% --- Perform a Standard (Two-Tailed) One-Sample T-test ---
% This is the most basic and stable function call.
fprintf('Performing a standard two-tailed t-test...\n');
t_results_2tailed = ttest(group_q_stats_data);

% --- Manually Calculate One-Tailed P-values for a Right-Tailed Test ---
% This is the robust way to get the result we want.
fprintf('Manually calculating right-tailed p-values...\n');

% Get the t-statistics and two-tailed p-values from the output object.
t_stats = t_results_2tailed.dat;
p_2tailed = t_results_2tailed.p;

% Create a new variable for one-tailed p-values.
p_1tailed = zeros(size(p_2tailed));

% For voxels with t > 0 (our direction of interest), p_1tailed is p_2tailed / 2.
positive_t_mask = t_stats > 0;
p_1tailed(positive_t_mask) = p_2tailed(positive_t_mask) ./ 2;

% For voxels with t <= 0, the result is not significant for a right-tailed test.
% The p-value is 1 - (p_2tailed / 2). These will be large and not survive.
negative_t_mask = t_stats <= 0;
p_1tailed(negative_t_mask) = 1 - p_2tailed(negative_t_mask) ./ 2;

% --- Create a New Statistic Image with the Correct One-Tailed P-values ---
% We copy the original object and overwrite the p-values.
t_results = t_results_2tailed;
t_results.p = p_1tailed;
t_results.dat_descrip = 'Group one-sample t-test (manually corrected to right-tailed)';

% --- Save the Unthresholded Group Maps ---
fprintf('Saving unthresholded group-level maps...\n');
group_t_map_filename = fullfile(output_dir, 'group_t_map_unthresholded.nii');
write(t_results, 'fname', group_t_map_filename, 'overwrite'); % Saves t-stats

p_map_to_save = t_results; % Use our new object as a template
p_map_to_save.dat = t_results.p; % Put the one-tailed p-values in the data field
p_map_to_save.dat_descrip = 'Group one-sample t-test (right-tailed p-values)';
group_p_map_filename = fullfile(output_dir, 'group_p_map_unthresholded.nii');
write(p_map_to_save, 'fname', group_p_map_filename, 'overwrite');

% --- Correct for Multiple Comparisons & Save Thresholded Map ---
fprintf('Applying FDR correction (q < 0.05) and saving thresholded map...\n');
t_results_fdr_corrected = threshold(t_results, .05, 'fdr');

group_t_map_fdr_filename = fullfile(output_dir, 'group_t_map_fdr_05.nii');
write(t_results_fdr_corrected, 'fname', group_t_map_fdr_filename, 'overwrite');

% --- Display Final Results ---
fprintf('Displaying FDR-corrected results. Close the figure to end the script.\n');
orthviews(t_results_fdr_corrected);

fprintf('\n--- Analysis Complete ---\n');
fprintf('You can find all output NIfTI files in: %s\n', output_dir);