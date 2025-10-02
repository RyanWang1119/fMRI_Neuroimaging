%% Load data
dat = fmri_data('data\WM_100307_res4d.nii');
time_series_data = dat.dat;

%% Misfit severity map
% From the 2008 paper.
% --- Setup ---
n_voxels_total = size(time_series_data, 1);
sw_bandwidth = [3, 5, 7];
sw_kernel_type = 'uniform'; % 'uniform' / 'gaussian'

% --- Data Cleaning ---
valid_voxels = ~any(isnan(time_series_data), 2) & var(time_series_data, [], 2) > 0;
data_valid = time_series_data(valid_voxels, :);
n_voxels_valid = size(data_valid, 1);
fprintf('Found %d valid voxels out of %d total voxels.\n', n_voxels_valid, n_voxels_total);


% SW statistic is calculated on mean-centered residuals, 
preprocessed_residuals =  detrend(data_valid')';

% --- Vectorized Calculation ---
sw_values_valid_voxels = zeros(n_voxels_valid, length(sw_bandwidth));

for i = 1:length(sw_bandwidth)
    w = sw_bandwidth(i);
    window_size = 2*w + 1;
    
    if strcmpi(sw_kernel_type, 'uniform')
        kernel = ones(1, window_size) / sqrt(window_size);
    elseif strcmpi(sw_kernel_type, 'gaussian')
        kernel = normpdf(-w:w, 0, w/2); 
        kernel = kernel / sqrt(sum(kernel.^2));
    end
    
    Yw_valid_voxels = conv2(preprocessed_residuals, kernel, 'same');
    sw_values_valid_voxels(:, i) = max(abs(Yw_valid_voxels), [], 2);
end

% Find the max statistic across bandwidths for each valid voxel
sw_map_valid = max(sw_values_valid_voxels, [], 2);

% --- Reconstruct Full Brain Image ---
% Create the final full-brain map.
sw_map_final = nan(n_voxels_total, 1);

% Place the results for the valid voxels back into the correct locations
sw_map_final(valid_voxels) = sw_map_valid;

sw_obj = dat;
sw_obj.dat = sw_map_final;
sw_obj.dat_descrip = 'Misfit Severity Map (SW Statistic)';

%% Plotting
orthviews(sw_obj);
title('Misfit Severity Map');