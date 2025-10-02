%% Load data
dat = fmri_data('data\WM_100307_res4d.nii');
time_series_data = dat.dat;

%% --- Spectral Analysis ---
TR = 2.0;  % might need to change
n_time_points = size(time_series_data, 2);
n_voxels = size(time_series_data, 1);

% Preprocessing
% Remove voxels with zero variance or NaN 
valid_voxels = ~any(isnan(time_series_data), 2) & var(time_series_data, [], 2) > 0;
fprintf('Found %d valid voxels out of %d\n', sum(valid_voxels), n_voxels);
% Detrend and normalize
preprocessed_residuals = detrend(time_series_data')';
preprocessed_residuals = preprocessed_residuals ./ std(preprocessed_residuals, [], 2); % normalize

% Apply windowing
win = hann(n_time_points)';
windowed_residuals = preprocessed_residuals .* win;
% Perform FFT
power_spectrum = abs(fft(windowed_residuals, [], 2)).^2 / sum(win.^2);
half_spectrum = power_spectrum(:, 1:floor(n_time_points/2)+1);

% Create frequency vector
fs = 1/TR;  % Sampling frequency
freqs = (0:floor(n_time_points/2)) * fs / n_time_points;

% Filter the frequency
min_freq = 0.01; 
max_freq = 0.25;  
freq_idx = find(freqs >= min_freq & freqs <= max_freq);

search_spectrum = half_spectrum(:, freq_idx);
search_freqs = freqs(freq_idx);

% Find dominant peaks
[peak_power, peak_idx_in_search] = max(search_spectrum, [], 2);
mean_power = mean(search_spectrum, 2);

% More sophisticated significance test
noise_floor = prctile(search_spectrum, 95, 2);  % 95th percentile as noise floor
is_significant = peak_power > 2 * noise_floor & valid_voxels;
fprintf('Found %d voxels with significant periodic activity\n', sum(is_significant));

%% Create result maps
periodicity_map_TRs = nan(n_voxels, 1);
frequency_map_Hz = nan(n_voxels, 1);
power_map = nan(n_voxels, 1);

% Only fill significant voxels
sig_idx = find(is_significant);
peak_freq_idx = peak_idx_in_search(is_significant);
dominant_freqs = search_freqs(peak_freq_idx);

frequency_map_Hz(sig_idx) = dominant_freqs;
periodicity_map_TRs(sig_idx) = 1 ./ dominant_freqs / TR;  % Convert to TRs
power_map(sig_idx) = peak_power(is_significant);

% Create fmri_data objects
freq_obj = dat;
freq_obj.dat = frequency_map_Hz;
freq_obj.dat_descrip = 'Dominant Frequency (Hz)';

period_obj = dat;
period_obj.dat = periodicity_map_TRs;
period_obj.dat_descrip = 'Dominant Period (TRs)';

power_obj = dat;
power_obj.dat = power_map;
power_obj.dat_descrip = 'Peak Spectral Power';

%% Plotting
% 1. Frequency Map
valid_freq_vals = frequency_map_Hz(~isnan(frequency_map_Hz));
freq_obj_thresh = threshold(freq_obj, [min(valid_freq_vals), max(valid_freq_vals)], 'raw-between');
orthviews(freq_obj_thresh);
title('Dominant Frequency Map (Hz)');
  
% 2. Period Map
valid_period_vals = periodicity_map_TRs(~isnan(periodicity_map_TRs));
period_obj_thresh = threshold(period_obj, [min(valid_period_vals), max(valid_period_vals)], 'raw-between');
orthviews(period_obj_thresh);
title('Dominant Period Map (TRs)');

% 3. Power Map
valid_power_vals = power_map(~isnan(power_map));
power_obj_thresh = threshold(power_obj, [prctile(valid_power_vals, 75), max(valid_power_vals)], 'raw-between');
orthviews(power_obj_thresh);
title('Peak Spectral Power');
