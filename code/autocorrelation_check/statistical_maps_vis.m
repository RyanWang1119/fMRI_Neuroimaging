%% --- Load Data ---
files_to_load = {
    'statistical_maps_WM_100307_res4d.mat', ...
    'statistical_maps_WM_108121_res4d.mat', ...
    'statistical_maps_WM_200917_res4d.mat'
    };

alpha = 0.05;

num_subjects = length(files_to_load);
results = struct(); 

for i = 1:num_subjects
    filename = files_to_load{i};   
    loaded_data = load(filename);
    subject_id = extractBetween(filename, 'statistical_maps_WM_', '_res4d.mat');
    
    results(i).subject_id = subject_id{1};
    results(i).lbq_obj = loaded_data.lbq_obj;
    results(i).adf_obj = loaded_data.adf_obj;
end

%% --- Plotting and Interpretation ---
% 1. Ljung-Box Map: Visualize voxels with significant autocorrelation
lbq_obj_thresh_1 = threshold(results(1).lbq_obj, [0, alpha], 'raw-between');
orthviews(lbq_obj_thresh_1);
title(sprintf('Voxels with Significant Autocorrelation (Ljung-Box p < %.2f)', alpha));

lbq_obj_thresh_2 = threshold(results(2).lbq_obj, [0, alpha], 'raw-between');
orthviews(lbq_obj_thresh_2);
title(sprintf('Voxels with Significant Autocorrelation (Ljung-Box p < %.2f)', alpha));

lbq_obj_thresh_3 = threshold(results(2).lbq_obj, [0, alpha], 'raw-between');
orthviews(lbq_obj_thresh_3);
title(sprintf('Voxels with Significant Autocorrelation (Ljung-Box p < %.2f)', alpha));

% 2. ADF Map: Visualize voxels that are likely stationary
adf_obj_thresh_1 = threshold(results(1).adf_obj, [0, alpha], 'raw-between');
orthviews(adf_obj_thresh_1);
title(sprintf('Likely Stationary Voxels (ADF p < %.2f)', alpha));


%% --- PLOT 2: P-Value Distribution Histograms ---
figure('Name', 'Ljung-Box P-value Distributions');
hold on;
for i = 1:num_subjects
    p_values = results(i).lbq_obj.dat(~isnan(results(i).lbq_obj.dat));
    histogram(p_values, 100, 'DisplayStyle', 'stairs', 'LineWidth', 1.5, ...
        'DisplayName', ['Subject ' results(i).subject_id]);
end
hold off;
title('Distribution of Ljung-Box P-values (Autocorrelation)');
xlabel('P-value');
ylabel('Number of Voxels');
legend;
grid on;
xline(alpha, 'r--', 'LineWidth', 2, 'Label', sprintf('alpha = %.2f', alpha));


%% --- PLOT 3: Ljung-Box Overlap Map  ---
custom_colors = [
    0 0 1;  % Blue for 1 subject
    0 1 0;  % Green for 2 subjects
    1 0 0;  % Red for 3 subjects
];

lbq_sig_maps = zeros(size(results(1).lbq_obj.dat, 1), num_subjects);
for i = 1:num_subjects
    lbq_sig_maps(:, i) = results(i).lbq_obj.dat < alpha & ~isnan(results(i).lbq_obj.dat);
end
lbq_overlap_data = sum(lbq_sig_maps, 2);
lbq_overlap_obj = results(1).lbq_obj; 
lbq_overlap_obj.dat = lbq_overlap_data;
lbq_overlap_obj.dat_descrip = 'Overlap map: Number of subjects with significant Ljung-Box test';

figure('Name', 'Ljung-Box Overlap Map');
orthviews(lbq_overlap_obj, 'clim', [1 3]);
colormap(custom_colors);
title('Ljung-Box: Voxel Overlap Across Subjects');

h_ax_leg = axes('Position', [0 0 1 1], 'Visible', 'off');
hold(h_ax_leg, 'on');

% build a legend manually.
p1 = plot(h_ax_leg, NaN, NaN, 's', 'MarkerFaceColor', custom_colors(1,:), 'MarkerEdgeColor', 'none', 'MarkerSize', 12);
p2 = plot(h_ax_leg, NaN, NaN, 's', 'MarkerFaceColor', custom_colors(2,:), 'MarkerEdgeColor', 'none', 'MarkerSize', 12);
p3 = plot(h_ax_leg, NaN, NaN, 's', 'MarkerFaceColor', custom_colors(3,:), 'MarkerEdgeColor', 'none', 'MarkerSize', 12);

lgd = legend(h_ax_leg, [p1, p2, p3], ...
       {'1 Subject', '2 Subjects', '3 Subjects'}, ...
       'Location', 'northeast', 'FontSize', 12);
       
title(lgd, 'Overlap Count');
hold(h_ax_leg, 'off');

%% --- PLOT 4: ADF Overlap Map  ---
adf_sig_maps = zeros(size(results(1).adf_obj.dat, 1), num_subjects);
for i = 1:num_subjects
    adf_sig_maps(:, i) = results(i).adf_obj.dat < alpha & ~isnan(results(i).adf_obj.dat);
end
adf_overlap_data = sum(adf_sig_maps, 2);
adf_overlap_obj = results(1).adf_obj; 
adf_overlap_obj.dat = adf_overlap_data;


figure('Name', 'ADF Overlap Map');
orthviews(adf_overlap_obj, 'clim', [1 3]);
colormap(custom_colors);
title('ADF: Voxel Overlap Across Subjects');

h_ax_leg2 = axes('Position', [0 0 1 1], 'Visible', 'off');
hold(h_ax_leg2, 'on');
p1_2 = plot(h_ax_leg2, NaN, NaN, 's', 'MarkerFaceColor', custom_colors(1,:), 'MarkerEdgeColor', 'none', 'MarkerSize', 12);
p2_2 = plot(h_ax_leg2, NaN, NaN, 's', 'MarkerFaceColor', custom_colors(2,:), 'MarkerEdgeColor', 'none', 'MarkerSize', 12);
p3_2 = plot(h_ax_leg2, NaN, NaN, 's', 'MarkerFaceColor', custom_colors(3,:), 'MarkerEdgeColor', 'none', 'MarkerSize', 12);
lgd2 = legend(h_ax_leg2, [p1_2, p2_2, p3_2], ...
       {'1 Subject', '2 Subjects', '3 Subjects'}, ...
       'Location', 'northeast', 'FontSize', 12);
title(lgd2, 'Overlap Count');
hold(h_ax_leg2, 'off');