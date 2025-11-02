%% Load data
dat = fmri_data('data\WM_200917_res4d.nii');
time_series_data = double(dat.dat);

%% --- Setup for Tests ---
n_voxels_total = size(time_series_data, 1);
n_time_points = size(time_series_data, 2);
alpha = 0.05; 
lags_for_lbq = min(20, floor(log(n_time_points))); 

%% --- Data Cleaning ---
% Find voxels with valid data (non-constant and no NaNs)
valid_voxels = ~any(isnan(time_series_data), 2) & var(time_series_data, [], 2) > 1e-6;
data_valid = time_series_data(valid_voxels, :);
n_voxels_valid = size(data_valid, 1);
fprintf('Found %d valid voxels out of %d total voxels.\n', n_voxels_valid, n_voxels_total);
preprocessed_residuals = detrend(data_valid')';

%% --- Voxel-wise Statistical Tests ---
if isempty(gcp('nocreate'))
    parpool; 
end

% Pre-allocate result vectors for valid voxels
lbq_p_values = nan(n_voxels_valid, 1);
adf_p_values = nan(n_voxels_valid, 1);

parfor v = 1:n_voxels_valid
    voxel_series = preprocessed_residuals(v, :);
    % 1. Ljung-Box Q-test
    [~, pL] = lbqtest(voxel_series, 'Lags', lags_for_lbq);
    lbq_p_values(v) = pL;
    % 2. Augmented Dickey-Fuller test
    [~, pA] = adftest(voxel_series, 'model', 'AR');
    adf_p_values(v) = pA;
end

%% --- Reconstruct Full Brain Images ---
% Create final full-brain maps, initializing with NaNs
lbq_map = nan(n_voxels_total, 1);
adf_map = nan(n_voxels_total, 1);

lbq_map(valid_voxels) = lbq_p_values;
adf_map(valid_voxels) = adf_p_values;

%% --- Create fmri_data Objects ---
% Ljung-Box p-value map
lbq_obj = dat;
lbq_obj.dat = lbq_map;
lbq_obj.dat_descrip = 'Ljung-Box Test p-values (Low p = Autocorrelated)';

% ADF p-value map
adf_obj = dat;
adf_obj.dat = adf_map;
adf_obj.dat_descrip = 'ADF Test p-values (Low p = Stationary)';

%% --- Plotting and Interpretation ---
% 1. Ljung-Box Map: Visualize voxels with significant autocorrelation
lbq_obj_thresh = threshold(lbq_obj, [0, alpha], 'raw-between');
orthviews(lbq_obj_thresh);
title(sprintf('Voxels with Significant Autocorrelation (Ljung-Box p < %.2f)', alpha));

% 2. ADF Map: Visualize voxels that are likely stationary
adf_obj_thresh = threshold(adf_obj, [0, alpha], 'raw-between');
orthviews(adf_obj_thresh);
title(sprintf('Likely Stationary Voxels (ADF p < %.2f)', alpha));
