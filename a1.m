% Step 1: Set Up Parameters and Select Voxel Data
% ----------------------------------------------------

% --- Define Experiment Parameters (You must set these) ---
TR = 2; % Repetition Time in seconds (e.g., 2s)
stim_onsets_sec = [10, 20, 30, 40, 50, 60, 70, 80]; % EXAMPLE onsets in seconds
fir_window_sec = 20; % How many seconds post-stimulus to model the error (e.g., 20s)

% --- Select Data ---
% Your residual data matrix
% Assuming dimensions are (Time, X, Y, Subject)
% Results = ...; 

% Extract data for one voxel for Subject 1
x = 1; y = 1; subj = 1;
voxel_residuals = Results(:, x, y, subj);

% --- Convert units from seconds to TRs (time points) ---
num_time_points = size(voxel_residuals, 1);
stim_onsets_tr = round(stim_onsets_sec / TR);
fir_window_tr = round(fir_window_sec / TR);


% Step 2: Create the FIR Design Matrix
% ----------------------------------------------------
% This matrix will have one column for each time point in the FIR window.
% Each row corresponds to a time point in the scan.
% This method creates a flexible basis set to estimate an arbitrary shape[cite: 1357].

fir_design_matrix = zeros(num_time_points, fir_window_tr);

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


% Step 3: Fit the Linear Model to the Residuals
% ----------------------------------------------------
% We are fitting the model: Residuals = FIR_Design_Matrix * Beta_Error + Noise
% The 'regress' function is useful as it provides confidence intervals.

% Add a constant term (intercept) to the design matrix
X = [ones(num_time_points, 1), fir_design_matrix];

[beta_error, beta_ci] = regress(voxel_residuals, X);

% The first coefficient is the intercept; the rest are for the FIR time points.
fir_coeffs = beta_error(2:end);
fir_ci = beta_ci(2:end, :);


% Step 4: Plot the Results
% ----------------------------------------------------
% Plot the FIR coefficients to visualize the shape of the systematic error.

time_axis_sec = (0:fir_window_tr-1) * TR;

figure;
hold on;

% Plot a horizontal line at zero for reference
plot([time_axis_sec(1), time_axis_sec(end)], [0, 0], 'k--'); 

% Plot the confidence interval as a shaded area
fill([time_axis_sec, fliplr(time_axis_sec)], [fir_ci(:,1)', fliplr(fir_ci(:,2)')], ...
     [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% Plot the estimated coefficients
plot(time_axis_sec, fir_coeffs, 'b-', 'LineWidth', 2);

hold off;

title('Estimated Shape of Systematic Error in Residuals');
xlabel('Time Since Stimulus Onset (s)');
ylabel('Magnitude of Unmodeled Response (beta)');
legend({'', '95% Confidence Interval', 'Estimated Error Shape'});
grid on;