load('data\task_residual\WMcHRFderiv.mat')

%% Mean of mean
% dim: [regions, time_points, subjects]
contrast_body   = squeeze(Results(:,:,:,5) - Results(:,:,:,1));
contrast_faces  = squeeze(Results(:,:,:,6) - Results(:,:,:,2));
contrast_places = squeeze(Results(:,:,:,7) - Results(:,:,:,3));
contrast_tools  = squeeze(Results(:,:,:,8) - Results(:,:,:,4));

% --- Average across all Subjects ---
mean_body   = mean(contrast_body, 3, 'omitnan');
mean_faces  = mean(contrast_faces, 3, 'omitnan');
mean_places = mean(contrast_places, 3, 'omitnan');
mean_tools  = mean(contrast_tools, 3, 'omitnan');

mean_body_roi   = mean(contrast_body, 1, 'omitnan');
mean_faces_roi  = mean(contrast_faces, 1, 'omitnan');
mean_places_roi = mean(contrast_places, 1, 'omitnan');
mean_tools_roi  = mean(contrast_tools, 1, 'omitnan');

% --- Average across all ROI ---
plot_body   = mean(mean_body, 1, 'omitnan');
plot_faces  = mean(mean_faces, 1, 'omitnan');
plot_places = mean(mean_places, 1, 'omitnan');
plot_tools  = mean(mean_tools, 1, 'omitnan');

figure;
hold on;
plot(plot_body, 'LineWidth', 2, 'DisplayName', 'Body (2bk-0bk)');
plot(plot_faces, 'LineWidth', 2, 'DisplayName', 'Faces (2bk-0bk)');
plot(plot_places, 'LineWidth', 2, 'DisplayName', 'Places (2bk-0bk)');
plot(plot_tools, 'LineWidth', 2, 'DisplayName', 'Tools (2bk-0bk)');
hold off;
title('Mean Residual (Across all Subjects and Regions)', 'FontSize', 20);
xlabel('Time Points');
ylabel('Mean Residual Amplitude');
legend('show', 'Location', 'best');
grid on;
yline(0, 'k--', 'Zero');

%% Mean ± Regional Standard Deviation
% Calculate Mean and SD across subjects ---
line_body = mean(mean_body_roi, 3, 'omitnan');
se_body  = std(mean_body_roi, 0, 3, 'omitnan')/sqrt(393);

line_faces = mean(mean_faces_roi, 3, 'omitnan');
se_faces  = std(mean_faces_roi, 0, 3, 'omitnan')/sqrt(393);

line_places = mean(mean_places_roi, 3, 'omitnan');
se_places  = std(mean_places_roi, 0, 3, 'omitnan')/sqrt(393);

line_tools = mean(mean_tools_roi, 3, 'omitnan');
se_tools  = std(mean_tools_roi, 0, 3, 'omitnan')/sqrt(393);


figure(2);
t = (1:size(line_body, 2))';
colors = get(gca, 'ColorOrder');

subplot(2, 2, 1);
plot_with_error(t, line_body, se_body, colors(1,:), 'Body');
title('Body');
grid on;
yline(0, 'k--');
xlabel('Time Point');
ylabel('Mean Residual');
subplot(2, 2, 2);
plot_with_error(t, line_faces, se_faces, colors(2,:), 'Faces');
title('Faces');
grid on;
yline(0, 'k--');
xlabel('Time Point');
ylabel('Mean Residual');
subplot(2, 2, 3);
plot_with_error(t, line_places, se_places, colors(3,:), 'Places');
title('Places');
grid on;
yline(0, 'k--');
xlabel('Time Point');
ylabel('Mean Residual');
subplot(2, 2, 4);
plot_with_error(t, line_tools, se_tools, colors(4,:), 'Tools');
title('Tools');
grid on;
yline(0, 'k--');
xlabel('Time Point');
ylabel('Mean Residual');

sgtitle('Mean Residuals of the Whole Brain (95% CI)', 'FontSize', 18, 'FontWeight', 'bold');


function plot_with_error(x, y_mean, y_error, color, display_name)
    fill_x = [x; flipud(x)];  
    fill_y = [y_mean' - 2 * y_error'; flipud(y_mean' + 2 * y_error')];
    
    fill(fill_x, fill_y, color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    hold on;
    plot(x, y_mean, 'Color', color, 'LineWidth', 2, 'DisplayName', display_name);
    hold off;
end

%% Correlation matrix
data_matrix = [plot_body', plot_faces', plot_places', plot_tools'];
contrast_names = {'Body', 'Faces', 'Places', 'Tools'};
R = corr(data_matrix);

disp('Trajectory Correlation Matrix (R):');
disp(R);

figure;
h = heatmap(contrast_names, contrast_names, R);
h.Title = 'Correlation of Grand Average Trajectories';
h.XLabel = 'Contrast';
h.YLabel = 'Contrast';
h.ColorLimits = [-1, 1];

%% Similarity tests on mean residuals (across subjects and regions): dTW, xcorrelation, Frechet
num_series = size(data_matrix, 2);
dtw_matrix = zeros(num_series, num_series);
peak_lags_matrix = zeros(num_series, num_series);
fre_matrix = zeros(num_series, num_series);

for i = 1:num_series
    for j = 1:num_series
        % DTW distance
        dtw_matrix(i, j) = dtw(data_matrix(:, i), data_matrix(:, j));
        % Cross-correlation
        [c, lags] = xcorr(data_matrix(:, i), data_matrix(:, j));
        [~, max_index] = max(c);
        peak_lags_matrix(i, j) = lags(max_index);
        % Frechet Distance
        fre_matrix(i, j) = frechetDistance(data_matrix(:, i), data_matrix(:, j));
    end
end

%% Dynamtic time wraps
figure;
h = heatmap(contrast_names, contrast_names, dtw_matrix);
h.Title = 'Pairwise DTW Distance Matrix';
h.XLabel = 'Time Series';
h.YLabel = 'Time Series';
h.ColorbarVisible = 'on';

%% Dynamtic time wraps
figure;
h = heatmap(contrast_names, contrast_names, peak_lags_matrix);
h.Title = 'Pairwise Lag Difference Matrix';
h.XLabel = 'Time Series';
h.YLabel = 'Time Series';
h.ColorbarVisible = 'on';

%% Fréchet Distance
figure;
h = heatmap(contrast_names, contrast_names, fre_matrix);
h.Title = 'Pairwise Fréchet Distance Matrix';
h.XLabel = 'Time Series';
h.YLabel = 'Time Series';
h.ColorbarVisible = 'on';

%% bar plot
n = length(contrast_names);
dtw_distances = [];
peak_lags = [];
fre_distances = [];
labels = {};
k = 1; 
for i = 1:n
    for j = (i+1):n 
        dtw_distances(k) = dtw_matrix(i, j);
        peak_lags(k) = dtw_matrix(i, j);
        fre_distances(k) = fre_matrix(i, j);
        labels{k} = [contrast_names{i}, ' - ', contrast_names{j}];
        
        k = k + 1;
    end
end

[sorted_dtw, sort_indices] = sort(dtw_distances, 'descend');
sorted_labels = labels(sort_indices);
sorted_peak_lags = peak_lags(sort_indices);
sorted_fre = fre_distances(sort_indices);

figure; 

ax1 = subplot(3, 1, 1);
bar(sorted_dtw);
title('DTW Distances');
ylabel('DTW Distance');
grid on;
ax1.XTick = 1:length(sorted_labels);
ax1.XTickLabel = [];

ax2 = subplot(3, 1, 2);
bar(sorted_peak_lags);
title('Peak Lags');
ylabel('Peak Lag Value');
grid on;
ax2.XTick = 1:length(sorted_labels);
ax2.XTickLabel = []; 

ax3 = subplot(3, 1, 3); 
bar(sorted_fre);
title('Frechet Distances');
ylabel('Frechet Distance');
grid on;
ax3.XTick = 1:length(sorted_labels);
ax3.XTickLabel = sorted_labels;
xtickangle(ax3, 45); 

%% Multivariate Distance Comparison (Regions as Dimensions)
data_cell_multi = {mean_body', mean_faces', mean_places', mean_tools'};
contrast_names = {'Body', 'Faces', 'Places', 'Tools'}; 
num_series = length(data_cell_multi);

fre_matrix_multi = zeros(num_series, num_series);
dtw_matrix_multi = zeros(num_series, num_series);

for i = 1:num_series
    for j = 1:num_series
        
        P = data_cell_multi{i}; 
        Q = data_cell_multi{j};

        fre_matrix_multi(i, j) = frechetDistance(P, Q);
        dtw_matrix_multi(i, j) = dtw(P, Q);
    end
end

%% Multivariate Fréchet Distance Heatmap
figure;
h = heatmap(contrast_names, contrast_names, fre_matrix_multi);
h.Title = 'Pairwise Multivariate Fréchet Distance (Regions as Dimensions)';
h.XLabel = 'Time Series';
h.YLabel = 'Time Series';
h.ColorbarVisible = 'on';

%% Multivariate DTW Distance Heatmap
figure;
h = heatmap(contrast_names, contrast_names, dtw_matrix_multi);
h.Title = 'Pairwise Multivariate DTW Distance (Regions as Dimensions)';
h.XLabel = 'Time Series';
h.YLabel = 'Time Series';
h.ColorbarVisible = 'on';

%% bar plot, Multivariate
n = length(contrast_names);
dtw_distances_multi = [];
peak_lags = [];
fre_distances_multi = [];
labels = {};
k = 1; 
for i = 1:n
    for j = (i+1):n 
        dtw_distances_multi(k) = dtw_matrix_multi(i, j);
        fre_distances_multi(k) = fre_matrix_multi(i, j);
        labels{k} = [contrast_names{i}, ' - ', contrast_names{j}];
        
        k = k + 1;
    end
end

[sorted_dtw_multi, sort_indices] = sort(dtw_distances_multi, 'descend');
sorted_labels = labels(sort_indices);
sorted_fre_multi = fre_distances_multi(sort_indices);

figure; 

ax1 = subplot(2, 1, 1);
bar(sorted_dtw);
title('DTW Distances');
ylabel('DTW Distance');
grid on;
ax1.XTick = 1:length(sorted_labels);
ax1.XTickLabel = [];

ax2 = subplot(2, 1, 2); 
bar(sorted_fre_multi);
title('Frechet Distances');
ylabel('Frechet Distance');
grid on;
ax2.XTick = 1:length(sorted_labels);
ax2.XTickLabel = sorted_labels;
xtickangle(ax2, 45); 

%%
Contrasts{1} = squeeze(Results(:,:,:,5) - Results(:,:,:,1)); % Body
Contrasts{2} = squeeze(Results(:,:,:,6) - Results(:,:,:,2)); % Faces
Contrasts{3} = squeeze(Results(:,:,:,7) - Results(:,:,:,3)); % Places
Contrasts{4} = squeeze(Results(:,:,:,8) - Results(:,:,:,4)); % Tools

titles = {'Body', 'Faces', 'Places', 'Tools'};
N = size(Contrasts{1}, 3);
t_vec = 1:size(Contrasts{1}, 2);

figure('Position', [100, 100, 1200, 800]);

for c = 1:4
    subplot(2, 2, c);
    hold on;

    mean_response = mean(Contrasts{c}, 3, 'omitnan');
    region_strength = sum(mean_response, 2, 'omitnan');
    
    [sorted_vals, sort_idx] = sort(region_strength, 'descend');
    
    top_idxs = sort_idx(5:10); mid_idxs = sort_idx(244:262); bot_idxs = sort_idx(478:489);  
       
    for i = 1:1
        r = top_idxs(i);
        data_roi = squeeze(Contrasts{c}(r, :, :));
        mu  = mean(data_roi, 2, 'omitnan'); sem = std(data_roi, 0, 2, 'omitnan') / sqrt(N);
        plot_with_error2(t_vec', mu, sem, [0.8 0.2 0.2]); 
    end
    
    for i = 1:1
        r = mid_idxs(i);
        data_roi = squeeze(Contrasts{c}(r, :, :));
        mu  = mean(data_roi, 2, 'omitnan'); sem = std(data_roi, 0, 2, 'omitnan') / sqrt(N);
        plot_with_error2(t_vec', mu, sem, [0.2 0.8 0.2]);
    end

    for i = 1:1
        r = bot_idxs(i);
        data_roi = squeeze(Contrasts{c}(r, :, :));
        mu  = mean(data_roi, 2, 'omitnan'); sem = std(data_roi, 0, 2, 'omitnan') / sqrt(N);
        plot_with_error2(t_vec', mu, sem, [0.2 0.2 0.8]);
    end
    
    yline(0, 'k--', 'LineWidth', 1.5);
    title(titles{c}, 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Time Point');
    ylabel('Residual Amplitude');
    if c == 1
        h1 = plot(nan, nan, 'Color', [0.8 0.2 0.2], 'LineWidth', 2);
        h2 = plot(nan, nan, 'Color', [0.2 0.8 0.2], 'LineWidth', 2);
        h3 = plot(nan, nan, 'Color', [0.2 0.2 0.8], 'LineWidth', 2);
        legend([h1 h2 h3], {'Top 3 "Underestimated"', 'Top 3 "Overestimated"'}, 'Location', 'best');
    end
    grid on;
    hold off;
end
sgtitle('Residual Trajectories of 3 Randomly Selected Regions (95% CI)', 'FontSize', 18);


% --- Helper Function ---
function plot_with_error2(x, y_mean, y_error, color)
    fill_x = [x; flipud(x)];
    fill_y = [y_mean - 2 * y_error; flipud(y_mean + 2 * y_error)];
    
    % Draw shaded area
    fill(fill_x, fill_y, color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    hold on;
    % Draw mean line
    plot(x, y_mean, 'Color', color, 'LineWidth', 1.5);
end

%% Spatial Heterogeneity
mean_body_regions   = mean(contrast_body, 3, 'omitnan');
mean_faces_regions  = mean(contrast_faces, 3, 'omitnan');
mean_places_regions = mean(contrast_places, 3, 'omitnan');
mean_tools_regions  = mean(contrast_tools, 3, 'omitnan');

% 4. Plotting
figure('Name', 'Spatial Heterogeneity of Residuals', 'Color', 'w');
t = 1:size(mean_body_regions, 2);

% -- Plot 1: Body --
subplot(2,2,1);
plot(t, mean_body_regions', 'Color', [0, 0.4470, 0.7410, 0.15]); 
hold on;
plot(t, mean(mean_body_regions, 1), 'k', 'LineWidth', 2); 
title('Body');
xlabel('Time Point'); ylabel('Residual Amplitude');
grid on; xlim([1 41]); yline(0, 'k--');

% -- Plot 2: Faces --
subplot(2,2,2);
plot(t, mean_faces_regions', 'Color', [0.8500, 0.3250, 0.0980, 0.15]);
hold on;
plot(t, mean(mean_faces_regions, 1), 'k', 'LineWidth', 2);
title('Faces');
xlabel('Time Point'); ylabel('Residual Amplitude');
grid on; xlim([1 41]); yline(0, 'k--');

% -- Plot 3: Places --
subplot(2,2,3);
plot(t, mean_places_regions', 'Color', [0.9290, 0.6940, 0.1250, 0.15]); 
hold on;
plot(t, mean(mean_places_regions, 1), 'k', 'LineWidth', 2);
title('Places');
xlabel('Time Point'); ylabel('Residual Amplitude');
grid on; xlim([1 41]); yline(0, 'k--');

% -- Plot 4: Tools --
subplot(2,2,4);
plot(t, mean_tools_regions', 'Color', [0.4940, 0.1840, 0.5560, 0.15]); 
hold on;
plot(t, mean(mean_tools_regions, 1), 'k', 'LineWidth', 2);
title('Tools');
xlabel('Time Point'); ylabel('Residual Amplitude');
grid on; xlim([1 41]); yline(0, 'k--');

sgtitle('Residual Trajectories of All 489 Regions', 'FontSize', 18);

%%
N = size(contrast_body, 3);
alpha = 0.05; % 95% CI
t_crit = tinv(1-alpha/2, N-1); % T-statistic for N-1 dof (approx 1.96 for large N)

% Helper function to calculate Mean and CI half-width
get_stats = @(x) deal( ...
    mean(squeeze(x), 2, 'omitnan'), ...           % Mean
    (std(squeeze(x), 0, 2, 'omitnan') / sqrt(N)) * t_crit ... % 95% CI
);

[m_body, ci_body]     = get_stats(mean_body_roi);
[m_faces, ci_faces]   = get_stats(mean_faces_roi);
[m_places, ci_places] = get_stats(mean_places_roi);
[m_tools, ci_tools]   = get_stats(mean_tools_roi);

% 5. Plotting
figure('Name', 'Global Mean Residuals Overlay', 'Color', 'w', 'Position', [100 100 900 600]);
clf; hold on;
t = (1:length(m_body))';

% Define Colors
colors = [0, 0.4470, 0.7410;    % Body (Blue)
          0.8500, 0.3250, 0.0980;  % Faces (Orange)
          0.9290, 0.6940, 0.1250;  % Places (Yellow)
          0.4940, 0.1840, 0.5560]; % Tools (Purple)

plot_ci_band(t, m_body, ci_body, colors(1,:));
plot_ci_band(t, m_faces, ci_faces, colors(2,:));
plot_ci_band(t, m_places, ci_places, colors(3,:));
plot_ci_band(t, m_tools, ci_tools, colors(4,:));

% -- Step B: Draw Mean Lines --
plot(t, m_body,   'Color', colors(1,:), 'LineWidth', 2, 'DisplayName', 'Body');
plot(t, m_faces,  'Color', colors(2,:), 'LineWidth', 2, 'DisplayName', 'Faces');
plot(t, m_places, 'Color', colors(3,:), 'LineWidth', 2, 'DisplayName', 'Places');
plot(t, m_tools,  'Color', colors(4,:), 'LineWidth', 2, 'DisplayName', 'Tools');

% -- Step C: Formatting --
hold off;
title('Mean Residuals of the Whole Brain (95% CI)', 'FontSize', 18);
xlabel('Time Point (0 to ~30s)');
ylabel('Mean Residual Amplitude');
yline(0, 'k--', 'Zero');
grid on;
legend('show', 'Location', 'best');
axis tight;

% --- Helper Function for Transparent CI Bands ---
function plot_ci_band(x, y_mean, y_ci, color)
    fill_x = [x; flipud(x)];
    fill_y = [y_mean - y_ci; flipud(y_mean + y_ci)];
    
    % Plot the patch with transparency
    h = fill(fill_x, fill_y, color, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    
    % Exclude this patch from the legend (so only lines show up)
    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
end