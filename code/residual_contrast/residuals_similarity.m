load('data\HRF_Resid_WMblock_LR.mat')

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
title('Average Residual (Mean across all subjects and regions)', 'FontSize', 14);
xlabel('Time Point (0 to ~30s)');
ylabel('Mean Residual Amplitude');
legend('show', 'Location', 'best');
grid on;
yline(0, 'k--', 'Zero');

%% Mean ± Regional Standard Deviation
% Calculate Mean and SD across Regions ---
line_body = mean(mean_body, 1, 'omitnan');
std_body  = std(mean_body, 0, 1, 'omitnan');

line_faces = mean(mean_faces, 1, 'omitnan');
std_faces  = std(mean_faces, 0, 1, 'omitnan');

line_places = mean(mean_places, 1, 'omitnan');
std_places  = std(mean_places, 0, 1, 'omitnan');

line_tools = mean(mean_tools, 1, 'omitnan');
std_tools  = std(mean_tools, 0, 1, 'omitnan');


figure(2);
t = (1:size(line_body, 2))';
colors = get(gca, 'ColorOrder');

subplot(2, 2, 1);
plot_with_error(t, line_body, std_body, colors(1,:), 'Body (2bk-0bk)');
title('Body (2bk-0bk)');
grid on;
yline(0, 'k--');
xlabel('Time Point');
ylabel('Mean Residual');
subplot(2, 2, 2);
plot_with_error(t, line_faces, std_faces, colors(2,:), 'Faces (2bk-0bk)');
title('Faces (2bk-0bk)');
grid on;
yline(0, 'k--');
xlabel('Time Point');
ylabel('Mean Residual');
subplot(2, 2, 3);
plot_with_error(t, line_places, std_places, colors(3,:), 'Places (2bk-0bk)');
title('Places (2bk-0bk)');
grid on;
yline(0, 'k--');
xlabel('Time Point');
ylabel('Mean Residual');
subplot(2, 2, 4);
plot_with_error(t, line_tools, std_tools, colors(4,:), 'Tools (2bk-0bk)');
title('Tools (2bk-0bk)');
grid on;
yline(0, 'k--');
xlabel('Time Point');
ylabel('Mean Residual');

sgtitle('Mean Residual (± 1 SD across Regions)', 'FontSize', 14, 'FontWeight', 'bold');


function plot_with_error(x, y_mean, y_error, color, display_name)
    fill_x = [x; flipud(x)];  
    fill_y = [y_mean' - y_error'; flipud(y_mean' + y_error')];
    
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