%% --- 0. Initial Setup ---
addpath(genpath('/users/rwang/CanLabCore'));
fprintf('--- Starting Voxel-wise ARMA Analysis ---\n');

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
% For each subject, concatenate the residual data from both LR and RL runs
file_LR = fullfile(base_dir, subject_id, 'LR', 'res4d.nii');
file_RL = fullfile(base_dir, subject_id, 'RL', 'res4d.nii');

if ~exist(file_LR, 'file') || ~exist(file_RL, 'file')
    fprintf('Data files for subject %s not found. Skipping.\n', subject_id);
    return;
end

fprintf('Loading data for subject %s...\n', subject_id);
data_LR = fmri_data(file_LR);
data_RL = fmri_data(file_RL);
combined_ts_data = [data_LR.dat; data_RL.dat];
[T, n_voxels] = size(combined_ts_data);
fprintf('Data loaded. Found %d voxels and %d total time points.\n', n_voxels, T);

%% --- 2. Voxel-wise ARMA Model Fitting ---
% --- Setup for parallel processing and model search ---
max_p = 3; 
max_q = 3; 
p_vals = 0:max_p;
q_vals = 0:max_q;
ar_order_map = zeros(1, n_voxels);
ma_order_map = zeros(1, n_voxels);

if isempty(gcp('nocreate'))
    parpool;
end
fprintf('Starting parallel ARMA model fitting for %d voxels...\n', n_voxels);
parfor v = 1:n_voxels
    voxel_ts = combined_ts_data(:, v);
    
    if var(voxel_ts) == 0
        continue;
    end
    
    best_bic = inf;
    best_p = 0;
    best_q = 0;
    
    % --- Grid search for the best (p,q) combination ---
    for p = p_vals
        for q = q_vals
            try
                Mdl = arima(p, 0, q);
                [~, ~, logL] = estimate(Mdl, voxel_ts, 'Display','off');
                num_params = p + q + 1; 
                [~, bic] = aicbic(logL, num_params, T);
                
                if bic < best_bic
                    best_bic = bic;
                    best_p = p;
                    best_q = q;
                end
                
            catch ME
            end
        end
    end
    ar_order_map(v) = best_p;
    ma_order_map(v) = best_q;
end

fprintf('Finished model fitting.\n');

%% --- 3. Save Results as NIFTI Brain Maps ---
% This corrected section first creates a proper 3D template object from the 
% 4D input data, then inserts the AR and MA order results before writing.

fprintf('Saving AR and MA order maps...\n');

% 1. Create a 3D template by selecting the first volume from the original data.
%    This ensures the new object's metadata correctly describes a 3D image.
template_3d_obj = select_volumes(data_LR, 1);

% --- Create and write the AR order map ---
ar_map_obj = template_3d_obj; % Copy the clean 3D template
ar_map_obj.dat = ar_order_map'; % Replace its data with the AR orders (note the transpose ')
ar_map_obj.fullpath = fullfile(output_dir, [subject_id '_ar_order_map.nii']);
write(ar_map_obj, 'overwrite');

% --- Create and write the MA order map ---
ma_map_obj = template_3d_obj; % Copy the clean 3D template
ma_map_obj.dat = ma_order_map'; % Replace its data with the MA orders (note the transpose ')
ma_map_obj.fullpath = fullfile(output_dir, [subject_id '_ma_order_map.nii']);
write(ma_map_obj, 'overwrite');

fprintf('--- Analysis for subject %s complete. Output maps saved in: %s ---\n', subject_id, output_dir);