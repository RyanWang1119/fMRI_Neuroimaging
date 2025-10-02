%% --- 0. Initial Setup ---
addpath(genpath('/users/rwang/CanLabCore'));
fprintf('--- Starting Step 1: Subject-Level Autocorrelation Mapping ---\n');

base_dir = '/dcs07/smart/data/HRF/bogdan_hrf/stats/fsl_noar_24vec_csf_spikes/results/WM';
output_dir = fullfile(base_dir, 'autocorrelation_analysis_output');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Basic check
all_items = dir(base_dir);
all_dirs = all_items([all_items.isdir]);
is_dot_dir = ismember({all_dirs.name}, {'.', '..'});
all_dirs = all_dirs(~is_dot_dir);

dir_names = {all_dirs.name};
is_numeric_6digit = ~cellfun('isempty', regexp(dir_names, '^\d{6}$'));
subject_dirs = all_dirs(is_numeric_6digit);

if isempty(subject_dirs)
    error('No 6-digit numeric subject directories found in %s. Please check the base_dir path.', base_dir);
end
subject_ids = {subject_dirs.name};
n_subjects = length(subject_ids);

task_id_str = getenv('SLURM_ARRAY_TASK_ID');
if isempty(task_id_str)
    error('This script is designed to run as a Slurm job array. SLURM_ARRAY_TASK_ID not found.');
end
i = str2double(task_id_str);

if i > n_subjects
    fprintf('Task ID %d is greater than the number of subjects (%d). Exiting.\n', i, n_subjects);
    return;
end

subject_id = subject_ids{i};
fprintf('Processing Subject %d/%d: %s\n', i, n_subjects, subject_id);

%% --- 1. Process LR and RL Runs for the Subject ---
% Define paths for both potential runs
file_LR = fullfile(base_dir, subject_id, 'LR', 'res4d.nii');
file_RL = fullfile(base_dir, subject_id, 'RL', 'res4d.nii');

qmaps_to_average = {}; % Cell array to hold the maps we generate
template_dat_obj = []; % To store header info for saving the final map

% --- Process LR Run ---
if exist(file_LR, 'file')
    fprintf('  Found LR run. Processing: %s\n', file_LR);
    [qmap, dat_obj] = calculate_qmap_for_run(file_LR);
    qmaps_to_average{end+1} = qmap;
    if isempty(template_dat_obj), template_dat_obj = dat_obj; end
else
    fprintf('  Warning: LR run not found for subject %s.\n', subject_id);
end

% --- Process RL Run ---
if exist(file_RL, 'file')
    fprintf('  Found RL run. Processing: %s\n', file_RL);
    [qmap, dat_obj] = calculate_qmap_for_run(file_RL);
    qmaps_to_average{end+1} = qmap;
    if isempty(template_dat_obj), template_dat_obj = dat_obj; end
else
    fprintf('  Warning: RL run not found for subject %s.\n', subject_id);
end

%% --- 2. Average Maps and Save Final Output ---

if isempty(qmaps_to_average)
    fprintf('Error: No valid res4d.nii files found for subject %s. Skipping.\n', subject_id);
    return;
end

all_qmaps = horzcat(qmaps_to_average{:});

% Average across the runs (the 2nd dimension), ignoring NaNs
final_q_map = nanmean(all_qmaps, 2);

fprintf('  Averaged %d map(s) for subject %s.\n', length(qmaps_to_average), subject_id);

% Create and Save the final fmri_data object
q_stat_obj = template_dat_obj; 
q_stat_obj.dat = final_q_map;
q_stat_obj.dat_descrip = 'Ljung-Box Q-statistic (averaged over LR/RL runs)';

output_q_filename = fullfile(output_dir, [subject_id '_Q_map.nii']);
write(q_stat_obj, 'fname', output_q_filename, 'overwrite');
fprintf('  Saved final averaged Q-statistic map to: %s\n', output_q_filename);

fprintf('\n--- Subject %s Complete ---\n', subject_id);


%% --- Helper Function ---
function [q_stat_map, dat_obj] = calculate_qmap_for_run(filepath)

    dat_obj = fmri_data(filepath);
    time_series_data = double(dat_obj.dat);

    % Setup for this run's test
    [n_voxels_total, n_time_points] = size(time_series_data);
    lags_for_lbq = min(10, floor(log(n_time_points)));

    % Data Cleaning
    valid_voxels_mask = var(time_series_data, 0, 2) > 1e-6 & ~any(isnan(time_series_data), 2);
    data_valid = time_series_data(valid_voxels_mask, :);
    n_voxels_valid = size(data_valid, 1);
    preprocessed_residuals = detrend(data_valid', 'constant')';

    % Voxel-wise Ljung-Box Test
    q_stats_vec = nan(n_voxels_valid, 1);
    for v = 1:n_voxels_valid
        [~, ~, stat, ~] = lbqtest(preprocessed_residuals(v, :), 'Lags', lags_for_lbq);
        q_stats_vec(v) = stat(end);
    end

    q_stat_map = nan(n_voxels_total, 1);
    q_stat_map(valid_voxels_mask) = q_stats_vec;
end