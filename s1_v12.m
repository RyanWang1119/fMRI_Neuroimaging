[n_time_points, dim_x, dim_y, n_subjects] = size(Results);

% --- Voxel for Residual Plot ---
voxel_to_plot = [1, 2]; 
subject_to_plot = 1;

%% --- Part 1: Plot Residual Time Course for a Single Voxel ---
% Extract the residual time series for the specified voxel and subject
residual_ts = squeeze(Results(:, voxel_to_plot(1), voxel_to_plot(2), subject_to_plot));

fig = figure('Name', 'Residual Time Course', 'NumberTitle', 'off', ...
             'Color', [0.95, 0.95, 0.95]); 
ax = axes('Parent', fig, ...
          'XColor', 'k', 'YColor', 'k', 'ZColor', 'k', ...
          'Color', [1, 1, 1], ... 
          'TickDir', 'in');
hold(ax, 'on');
% Plot the residual time course
plot(ax, residual_ts, 'b-');
plot(ax, zeros(1, n_time_points), 'r--');
title(ax, sprintf('Residual Time Series for Voxel [%d, %d], Subject %d', ...
                  voxel_to_plot(1), voxel_to_plot(2), subject_to_plot), ...
      'Color', 'k', 'FontSize', 14, 'FontWeight', 'bold');
xlabel(ax, 'Time (seconds)', 'Color', 'k');
ylabel(ax, 'Residual Value', 'Color', 'k');
legend(ax, 'Residuals', 'Zero Line');
grid(ax, 'on');

fig = figure('Name', 'Residual Histogram');
histogram(residual_ts, 50);
title('Histogram of residual\_ts');
xlabel('Value');
ylabel('Frequency');
grid on;

%% --- Part 2: ACF and PACF Plots ---
figure;
subplot(2,1,1);
autocorr(residual_ts);
title('Autocorrelation Function (ACF)');

subplot(2,1,2);
parcorr(residual_ts);
title('Partial Autocorrelation Function (PACF)');

%% Spectral Analysis Using Fourier Transformation
TR = 2; % Repetition Time in seconds. Need to change
subject_index = 1; % Analyzing the first subject

% Map 1: Stores the period in units of TRs 
periodicity_map_TRs = nan(dim_x, dim_y);
% Map 2: Stores the frequency in Hertz (Hz)
frequency_map_Hz = nan(dim_x, dim_y);

for x = 1:dim_x
    for y = 1:dim_y
        voxel_residuals = squeeze(Results(:, x, y, subject_index));
        
        preprocessed_residuals = detrend(voxel_residuals, 'constant'); 
        N = n_time_points;
        
        %  Apply a Hanning window to reduce spectral leakage.
        win = hann(N);
        windowed_residuals = preprocessed_residuals .* win;
        
        power_spectrum = abs(fft(windowed_residuals)).^2 / sum(win.^2);
        half_spectrum = power_spectrum(1:floor(N/2)+1);
        
        % Find Dominant Peak and Check for Significance 
        search_spectrum = half_spectrum(2:end);
        [peak_power, peak_idx_in_search] = max(search_spectrum); 
        mean_power = mean(search_spectrum);
        
        if peak_power > 3 * mean_power
            peak_idx = peak_idx_in_search + 1;
                       
            % Calculate the dominant period in TRs 
            dominant_period_in_TRs = N / (peak_idx - 1);
            % Calculate the dominant frequency in Hz
            frequency_in_Hz = (peak_idx - 1) / (N * TR);
            
            periodicity_map_TRs(x, y) = dominant_period_in_TRs;
            frequency_map_Hz(x, y) = frequency_in_Hz;
            
        end
    end
end


% Plot 1: Dominant Periodicity Map (in TRs)
figure('Name', 'Dominant Periodicity Map');
subplot(1, 2, 1);
imagesc(periodicity_map_TRs);
colorbar;
title({'Dominant Period (Seasonality)', sprintf('Subject %d (in TRs)', subject_index)});
xlabel('Y dimension');
ylabel('X dimension');


% Plot 2: Dominant Frequency Map (in Hz)
subplot(1, 2, 2);
imagesc(frequency_map_Hz);
colorbar;
title({'Dominant Frequency', sprintf('Subject %d (in Hz)', subject_index)});
xlabel('Y dimension');
ylabel('X dimension');

all_values = periodicity_map_TRs(:);
figure;
h = histogram(all_values);


%% Misfit severity map
% From the 2008 paper.
sw_bandwidth = [3, 5, 7];
sw_kernel_type = 'uniform'; % 'uniform' / 'gaussian'

sw_map = zeros(dim_x, dim_y);

for x = 1:dim_x
    for y = 1:dim_y
        voxel_residuals = squeeze(Results(:, x, y, subject_index));     
        sw_map(x, y) = calculateSwStatistic(voxel_residuals, sw_bandwidth, sw_kernel_type);
    end
end

function sw_max = calculateSwStatistic(residuals, bandwidths, kernel_type)
   
    residuals = detrend(residuals, 'constant'); 
    sw_values = zeros(length(bandwidths), 1); 
    for i = 1:length(bandwidths)
        w = bandwidths(i);
        window_size = 2*w + 1;

        if strcmpi(kernel_type, 'uniform')
            % Uniform kernel, normalized so that sum(K.^2) = 1
            kernel = ones(1, window_size) / sqrt(window_size);
        elseif strcmpi(kernel_type, 'gaussian')
            % Gaussian kernel, normalized so that sum(K.^2) = 1
            kernel = normpdf(-w:w, 0, w/2); 
            kernel = kernel / sqrt(sum(kernel.^2));
        end

        % Convolve to get the moving average of residuals
        Yw = conv(residuals, kernel, 'same');
        sw_values(i) = max(abs(Yw));
    end
    sw_max = max(sw_values);
end

figure('Name', 'Diagnostic Maps for fMRI Residuals', 'NumberTitle', 'off', 'Position', [100 100 1500 400]);
imagesc(sw_map);
colorbar;
title({'Misfit Severity Map', '(Sw Statistic)'});
xlabel('Y Dimension');
ylabel('X Dimension');
axis square;


%% Group-Level Model Misfit Analysis 
[n_time_points, dim_x, dim_y, num_subjects] = size(Results);

% Analysis Parameters
sw_bandwidths = [3, 5, 7]; 
sw_kernel_type = 'uniform';
mc_permutations = 1000;

% --- Start the Parallel Pool ---
if isempty(gcp('nocreate'))
    parpool;
end

% --- Step 1: Calculate p-values for all voxels and subjects ---
fprintf('Starting parallel p-value calculation for %d subjects...\n', num_subjects);
all_p_values = zeros(dim_x, dim_y, num_subjects);

for s = 1:num_subjects
    fprintf('Processing Subject %d...\n', s);
    
    current_subject_p_values = zeros(dim_x, dim_y);
    subject_data = Results(:, :, :, s);
    
    parfor x = 1:dim_x
       
        y_p_values = zeros(1, dim_y); 
        for y = 1:dim_y
            voxel_residuals = squeeze(subject_data(:, x, y));
            
            if var(voxel_residuals) == 0
                y_p_values(y) = 1;
                continue;
            end
            
            y_p_values(y) = calculatePValueMC(voxel_residuals, sw_bandwidths, sw_kernel_type, mc_permutations);
        end
        
        current_subject_p_values(x, :) = y_p_values;
    end
    
    all_p_values(:, :, s) = current_subject_p_values;
end

% --- Helper Functions ---
function p_val = calculatePValueMC(residuals, bandwidths, kernel_type, num_permutations)
    observed_sw = calculateSwStatistic(residuals, bandwidths, kernel_type);
    null_sw_distribution = zeros(num_permutations, 1);
    for i = 1:num_permutations
        permuted_residuals = residuals .* sign(randn(size(residuals)));
        null_sw_distribution(i) = calculateSwStatistic(permuted_residuals, bandwidths, kernel_type);
    end
    p_val = (sum(null_sw_distribution >= observed_sw) + 1) / (num_permutations + 1);
end

%% --- Plot ---
load("all_p_values.mat")
[dim_t, dim_x, dim_y, num_subjects] = size(Results);
sw_bandwidths = [3, 5, 7]; 
sw_kernel_type = 'uniform';
mc_permutations = 1000;
fprintf('Combining p-values across subjects...\n');
all_p_values(all_p_values == 0) = 1e-12; 
log_p_values = log(all_p_values);
Q_map = -2 * sum(log_p_values, 3);

% --- Step 3: Determine Significance and Plot ---
df = 2 * num_subjects;
p_thresh = 0.05;
q_critical = chi2inv(1 - p_thresh, df);
Q_map_thresholded = Q_map;
Q_map_thresholded(Q_map < q_critical) = 0;

% Plotting 
figure('Name', 'Group-Level Model Misfit Analysis', 'NumberTitle', 'off', 'Position', [100 100 1200 500]);
subplot(1, 2, 1);
imagesc(Q_map); colorbar;
title({'Group-Level Misfit Map (Q Statistic)', sprintf('Chi-Square(%d df)', df)});
xlabel('Y Dimension'); ylabel('X Dimension'); axis square;
subplot(1, 2, 2);
imagesc(Q_map_thresholded); colorbar;
title({'Significant Misfit', sprintf('p < %.2f (Q > %.2f)', p_thresh, q_critical)});
xlabel('Y Dimension'); ylabel('X Dimension'); axis square;



