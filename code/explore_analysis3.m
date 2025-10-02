dat = fmri_data('data\WM_100307_res4d.nii');

%% --- Spectral Analysis ---
TR = 2;  % need to change
time_series_data = dat.dat;
n_time_points = size(time_series_data, 2);
n_voxels = size(time_series_data, 1);

% detrend data
preprocessed_residuals = detrend(time_series_data')';

%  Apply a Hanning window to reduce spectral leakage.
win = hann(n_time_points)';
windowed_residuals = preprocessed_residuals .* win;

% Perform Fourier Transform to compute the power spectrum
power_spectrum = abs(fft(windowed_residuals, [], 2)).^2 / sum(win.^2);
half_spectrum = power_spectrum(:, 1:floor(n_time_points/2)+1);
search_spectrum = half_spectrum(:, 2:end);

% Find Dominant Peak and Check for Significance 
[peak_power, peak_idx_in_search] = max(search_spectrum, [], 2);
mean_power = mean(search_spectrum, 2);
is_significant = peak_power > 3 * mean_power;
periodicity_map_TRs = nan(n_voxels, 1);
frequency_map_Hz = nan(n_voxels, 1);

peak_idx = peak_idx_in_search(is_significant) + 1;
periodicity_map_TRs(is_significant) = n_time_points ./ (peak_idx - 1);
frequency_map_Hz(is_significant) = (peak_idx - 1) / (n_time_points * TR);

%% 3. Create fmri_data Objects for Results
freq_obj = dat;  
freq_obj.dat = frequency_map_Hz;  
freq_obj.dat_descrip = 'Dominant Frequency (Hz)';

period_obj = dat;  
period_obj.dat = periodicity_map_TRs; 
period_obj.dat_descrip = 'Dominant Period (TRs)';

%% 4. Visualize the Frequency Map 
min_freq = min(freq_obj.dat(freq_obj.dat > 0 & ~isnan(freq_obj.dat)));
max_freq = max(freq_obj.dat(~isnan(freq_obj.dat)));

fprintf('Applying threshold to create a new display object...\n');
freq_obj_thresh = threshold(freq_obj, [min_freq max_freq], 'raw-between'); 

orthviews(freq_obj_thresh)
colorbar;
