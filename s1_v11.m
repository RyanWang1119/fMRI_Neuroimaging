[n_time_points, dim_x, dim_y, n_subjects] = size(Results);

% --- Voxel for Residual Plot ---
voxel_to_plot = [1, 1]; 
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
%  ACF and PACF plots having significant spikes at lags 5, 10, 15 
%  and a marginally significant spike at lag 2, this suggests a seasonal pattern with period 5.

%% --- Part 3: Stationarity and White Noise Check ---
% Augmented Dickey-Fuller test (null: non-stationary)
[h, pValue, stat, cValue] = adftest(residual_ts);
if h == 1
    fprintf('ADF Test: Series is stationary (p-value: %.4f)\n', pValue);
else
    fprintf('ADF Test: Series is non-stationary (p-value: %.4f)\n', pValue);
end

% KPSS test (null: stationary)
[h, pValue] = kpsstest(residual_ts);
if h == 0
    fprintf('KPSS Test: Series is stationary (p-value: %.4f)\n', pValue);
else
    fprintf('KPSS Test: Series is non-stationary (p-value: %.4f)\n', pValue);
end

[h_lb, p_lb] = lbqtest(residual_ts, 'Lags', 1:10);
disp('Ljung-Box Test:');
fprintf('h = %d (1 = not white noise), p-value = %.4f\n', h_lb, p_lb);
% Small p-values for the initial lags. 
% The series is stationary and not a white noise.

%% --- Part 4: Modeling ---

% --- Model Estimation and Comparison ---
% Both ACF and PACF tailing off at the seasonal lags, probably a SARMA(1,1) model.

% Model 1: SARIMA(0,0,0)(1,0,0) with S=5
Mdl1 = arima( 'SARLags',5, 'Seasonality',5);
[EstMdl1, ~, logL1] = estimate(Mdl1, residual_ts, 'Display','off');
[aic1, bic1] = aicbic(logL1, 2, length(residual_ts)); 
% Model 2: SARIMA(0,0,0)(1,0,0) with S=5
Mdl2 = arima( 'SARLags',5, 'Seasonality',5);
[EstMdl2, ~, logL2] = estimate(Mdl2, residual_ts, 'Display','off');
[aic2, bic2] = aicbic(logL2, 2, length(residual_ts)); 
% Model 3: SARIMA(0,0,0)(1,0,1) with S=5
Mdl3 = arima('SARLags',5, 'SMALags',5, 'Seasonality',5);
[EstMdl3, ~, logL3] = estimate(Mdl3, residual_ts, 'Display','off');
[aic3, bic3] = aicbic(logL3, 2, length(residual_ts));

% results
fprintf('--- Model Comparison ---\n');
fprintf('Model                     |      AIC      |      BIC\n');
fprintf('----------------------------------------------------------\n');
fprintf('SARIMA(1,0,0)(1,0,0)5     |   %8.4f    |   %8.4f\n', aic1, bic1);
fprintf('SARIMA(0,0,1)(1,0,0)5     |   %8.4f    |   %8.4f\n', aic2, bic2);
fprintf('SARIMA(0,0,0)(1,0,1)5     |   %8.4f    |   %8.4f\n', aic3, bic3);


%% --- Step 5: Model Diagnostic Checking ---
EstMdl = estimate(Mdl3, residual_ts);
res = infer(EstMdl, residual_ts);

% 1. Plot the residuals
figure;
subplot(2,2,1);
plot(res);
title('Model Residuals');
grid on;

% 2. ACF of residuals
subplot(2,2,2);
autocorr(res);
title('ACF of Residuals');

% 3. PACF of residuals
subplot(2,2,3);
parcorr(res);
title('PACF of Residuals');

% 5. Ljung-Box Q-test for residual autocorrelation
% H0: The residuals are not autocorrelated 
[h_lbq, pValue_lbq] = lbqtest(res);
fprintf('--- Ljung-Box Test on Model Residuals ---\n');
fprintf('Ljung-Box Test Statistic (h): %d\n', h_lbq);
fprintf('p-value: %f\n', pValue_lbq);
if h_lbq == 0
    fprintf('Result: No significant residual autocorrelation.\n\n');
else
    fprintf('Result: Significant residual autocorrelation exists.\n\n');
end


%% --- Step 6: Non-seasonal part ---
% Model 4: SARIMA(1,0,0)(1,0,1) with S=5
Mdl4 = arima("ARLags",2, 'SARLags',5, 'SMALags',5, 'Seasonality',5);
[EstMdl4, ~, logL4] = estimate(Mdl4, residual_ts, 'Display','off');
[aic4, bic4] = aicbic(logL4, 2, length(residual_ts));

res = infer(EstMdl4, residual_ts);


% 1. Plot the residuals
figure;
subplot(2,2,1);
plot(res);
title('Model Residuals');
grid on;

% 2. ACF of residuals
subplot(2,2,2);
autocorr(res);
title('ACF of Residuals');

% 3. PACF of residuals
subplot(2,2,3);
parcorr(res);
title('PACF of Residuals');

% 5. Ljung-Box Q-test for residual autocorrelation
% H0: The residuals are not autocorrelated 
[h_lbq, pValue_lbq] = lbqtest(res);
fprintf('--- Ljung-Box Test on Model Residuals ---\n');
fprintf('Ljung-Box Test Statistic (h): %d\n', h_lbq);
fprintf('p-value: %f\n', pValue_lbq);
if h_lbq == 0
    fprintf('Result: No significant residual autocorrelation.\n\n');
else
    fprintf('Result: Significant residual autocorrelation exists.\n\n');
end

%% --- Step 7: Detecting systematic mis-modeling ---
% From the 2007 paper
% Use FIR to model the shape of the error.
% makes minimal assumptions about the shape of the mis-modeling
% Check if the error is systematically correlated with the experimental task.

% Set parameters, need to change
TR = 2; 
stim_onsets_sec = [30, 90, 150, 210, 270, 330, 390, 450]; 
fir_window_sec = 20; 
x = 1; y = 1; subj = 1;
voxel_residuals = Results(:, x, y, subj);

num_time_points = size(voxel_residuals, 1);
stim_onsets_tr = round(stim_onsets_sec / TR);
fir_window_tr = round(fir_window_sec / TR);

fir_design_matrix = zeros(num_time_points, fir_window_tr);

% Create the FIR design matrix
for i = 1:length(stim_onsets_tr)
    onset = stim_onsets_tr(i);
    
    for j = 1:fir_window_tr
        % The time point for the current lag (j) after the current onset (i)
        event_time_point = onset + j - 1;
        
        % Ensure the time point is within the scan duration
        if event_time_point <= num_time_points
            % The j-th column of the matrix corresponds to the j-th TR after an event.
            % Place a 1 at the corresponding time point.
            fir_design_matrix(event_time_point, j) = 1;
        end
    end
end

% Fit the Linear Model to the Residual
% Add a constant term (intercept) to the design matrix
X = [ones(num_time_points, 1), fir_design_matrix];
% Perform the linear regression, 
% the traditional least-square solution is very sensitive to noise, though (2007).
[beta_error, beta_ci] = regress(voxel_residuals, X);
fir_coeffs = beta_error(2:end);
fir_ci = beta_ci(2:end, :);


time_axis_sec = (0:fir_window_tr-1) * TR;

figure;
hold on;
plot([time_axis_sec(1), time_axis_sec(end)], [0, 0], 'k--', 'LineWidth', 1); 
fill([time_axis_sec, fliplr(time_axis_sec)], [fir_ci(:,1)', fliplr(fir_ci(:,2)')], ...
     [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(time_axis_sec, fir_coeffs, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'MarkerSize', 4);
hold off;

title('Estimated Shape of Systematic Error in Residuals');
xlabel('Time Since Stimulus Onset (s)');
ylabel('Magnitude of Unmodeled Response (Beta)');
legend({'', '95% Confidence Interval', 'Estimated Error Shape'}, 'Location', 'best');
grid on;
box on;

% The significant coefficients identify the locations where the model HRF and the actual HRF differ.
% The true BOLD signal was, on average, stronger than the model predicted

