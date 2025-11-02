%% Load results and atlas
load('data\WMblock_parallel_results.mat')  
addpath(genpath('C:\Users\Ryan_\CanlabCore'))
addpath('C:\Users\Ryan_\Neuroimaging_Pattern_Masks\Atlases_and_parcellations\2018_Wager_combined_atlas')
atlas_obj = load_atlas('canlab2018');
region_labels = atlas_obj.labels;  

time_points = 41;
contrasts = 6;
bsamp = 20;
regions = 489;
TR = 0.72;
time_vec = (0:time_points-1) * TR;  

contrast_names = {'body_vs_faces', 'body_vs_places', 'body_vs_tools', ...
                  'faces_vs_places', 'faces_vs_tools', 'places_vs_tools'};

%% 1. Compute Beta Statistics Across Bootstraps
% beta: feature importance
% Dimension 1 (size 41): To store the results for each time point.
% Dimension 2 (size 6): To store the results for each contrast.
% Dimension 3 (size 489): To store the 489 region weights that the SVM produced for that specific time point and contrast.

Beta_se = Beta_std / sqrt(bsamp);
% Consistency score 
Beta_consistency = abs(Beta_mean) ./ (Beta_se + eps);  % 41 x 6 x 489

%% 2. Identify Most Important Regions for Each Contrast

% Method 1: Average absolute weight across all time points
Beta_importance_overall = squeeze(mean(abs(Beta_mean), 1));  % 6 x 489

% Method 2: Peak absolute weight across time
Beta_importance_peak = squeeze(max(abs(Beta_mean), [], 1));  % 6 x 489

% Method 3: Average consistency score across time
Beta_consistency_overall = squeeze(mean(Beta_consistency, 1));  % 6 x 489

%% 3. Rank Top Regions for Each Contrast
top_n = 20;  

fprintf('\n=== TOP %d REGIONS BY CONTRAST ===\n', top_n);

Top_Regions = struct();

for c = 1:contrasts
    fprintf('\n--- %s ---\n', contrast_names{c});
    
    % Sort by overall importance (average absolute weight)
    [sorted_importance, sort_idx] = sort(Beta_importance_overall(c, :), 'descend');
    
    Top_Regions.(contrast_names{c}).indices = sort_idx(1:top_n);
    Top_Regions.(contrast_names{c}).importance = sorted_importance(1:top_n);
    Top_Regions.(contrast_names{c}).labels = region_labels(sort_idx(1:top_n));
    Top_Regions.(contrast_names{c}).mean_weights = Beta_mean(:, c, sort_idx(1:top_n));
    
   fprintf('%-6s %-15s %s\n', 'Rank', 'Region', 'Importance');
    fprintf('%-6s %-15s %s\n', '----', '------', '----------');
    for i = 1:top_n
        region_name = region_labels{sort_idx(i)};
        fprintf('%-6d %-15s %.4f\n', i, region_name, sorted_importance(i));
    end
end

%% 4. Visualize Region Importance Heatmaps
figure('Position', [100, 100, 1400, 800]);
sgtitle('Beta Weights Across Time and Regions', 'FontSize', 14, 'FontWeight', 'bold');

for c = 1:contrasts
    subplot(2, 3, c)
    
    % Get top regions for this contrast
    top_regions = Top_Regions.(contrast_names{c}).indices(1:top_n);
    
    % Extract beta weights for top regions only
    beta_subset = squeeze(Beta_mean(:, c, top_regions))';  % top_n x 41
    
    imagesc(time_vec, 1:top_n, beta_subset)
    colorbar
    
    % Create red-blue colormap (blue=negative, white=zero, red=positive)
    n = 256;
    redblue_cmap = [linspace(0, 1, n/2)', linspace(0, 1, n/2)', ones(n/2, 1); ...
                    ones(n/2, 1), linspace(1, 0, n/2)', linspace(1, 0, n/2)'];
    colormap(redblue_cmap)
    
    caxis([-max(abs(beta_subset(:))), max(abs(beta_subset(:)))])  % Symmetric color scale
    
    xlabel('Time (seconds)', 'FontSize', 10)
    ylabel('Top Regions (ranked)', 'FontSize', 10)
    title(strrep(contrast_names{c}, '_', ' '), 'FontSize', 11)
    
    % Add vertical line at typical HRF peak (~6 seconds)
    hold on
    plot([6 6], ylim, 'k--', 'LineWidth', 1.5)
    hold off
end

saveas(gcf, 'Region_Importance_Heatmap.png')

%% 5. Temporal Dynamics of Top Regions

figure('Position', [100, 100, 1400, 1000]);
sgtitle('Temporal Dynamics of Top 5 Regions per Contrast', 'FontSize', 14, 'FontWeight', 'bold');

for c = 1:contrasts
    subplot(3, 2, c)
    
    top_5_regions = Top_Regions.(contrast_names{c}).indices(1:5);
    top_5_labels = Top_Regions.(contrast_names{c}).labels(1:5);
    
    beta_top5 = squeeze(Beta_mean(:, c, top_5_regions));  % 41 x 5
    
    plot(time_vec, beta_top5, 'LineWidth', 2)
    xlabel('Time (seconds)', 'FontSize', 10)
    ylabel('Beta Weight', 'FontSize', 10)
    title(strrep(contrast_names{c}, '_', ' '), 'FontSize', 11)
    legend(top_5_labels, 'Location', 'best', 'FontSize', 8, 'Interpreter', 'none')
    grid on
    
    % Add horizontal line at zero
    hold on
    plot(xlim, [0 0], 'k--', 'LineWidth', 1)
    hold off
end

saveas(gcf, 'Top_Regions_Temporal_Dynamics.png')

%% 6. Region Overlap Analysis Across Contrasts

% Find regions that are important across multiple contrasts
threshold_percentile = 90;  % Top 10% of regions

Important_Regions_Per_Contrast = false(contrasts, regions);

for c = 1:contrasts
    thresh = prctile(Beta_importance_overall(c, :), threshold_percentile);
    Important_Regions_Per_Contrast(c, :) = Beta_importance_overall(c, :) > thresh;
end

% Count how many contrasts each region is important for
Region_Consistency = sum(Important_Regions_Per_Contrast, 1);  % 1 x 489

% Find regions important across multiple contrasts
multi_contrast_regions = find(Region_Consistency >= 3);  % Important in 3+ contrasts

fprintf('\n=== REGIONS IMPORTANT ACROSS MULTIPLE CONTRASTS ===\n');
fprintf('Regions important in 3 or more contrasts:\n\n');

[~, sort_idx] = sort(Region_Consistency(multi_contrast_regions), 'descend');
multi_contrast_regions_sorted = multi_contrast_regions(sort_idx);

for i = 1:length(multi_contrast_regions_sorted)
    reg_idx = multi_contrast_regions_sorted(i);
    fprintf('%s (important in %d contrasts)\n', ...
        region_labels{reg_idx}, Region_Consistency(reg_idx));
end

%% 7. Network-Level Analysis (if atlas has network labels)

% Check if atlas has network information
if isfield(atlas_obj, 'probability_maps')
    fprintf('\n=== NETWORK-LEVEL SUMMARY ===\n');
    
    % Extract network labels (simplified - adapt based on your atlas structure)
    % This assumes region labels contain network information
    % You may need to customize this based on your specific atlas
    
    for c = 1:contrasts
        fprintf('\n%s:\n', contrast_names{c});
        top_20_labels = Top_Regions.(contrast_names{c}).labels;
        
        % Count network representations in top regions
        network_counts = struct();
        for i = 1:length(top_20_labels)
            % Extract network name (assuming format like "Network_Region")
            parts = strsplit(top_20_labels{i}, '_');
            if length(parts) > 1
                network = parts{1};
                if isfield(network_counts, network)
                    network_counts.(network) = network_counts.(network) + 1;
                else
                    network_counts.(network) = 1;
                end
            end
        end
        
        % Display network summary
        network_names = fieldnames(network_counts);
        for j = 1:length(network_names)
            fprintf('  %s: %d regions\n', network_names{j}, network_counts.(network_names{j}));
        end
    end
end

%% 8. Statistical Significance Testing for Region Weights

% Test which regions have consistently non-zero weights across bootstraps
Alpha = 0.05;
Z_threshold = norminv(1 - Alpha/2);  % Two-tailed test with Bonferroni correction

Significant_Regions = struct();

for c = 1:contrasts
    % For each time point and region, test if significantly different from zero
    Z_scores = Beta_mean(:, c, :) ./ (Beta_se(:, c, :) + eps);  % 41 x 1 x 489
    
    % A region is significant if it's significant at ANY time point
    sig_any_time = any(abs(squeeze(Z_scores)) > Z_threshold, 1);  % 1 x 489
    
    Significant_Regions.(contrast_names{c}).indices = find(sig_any_time);
    Significant_Regions.(contrast_names{c}).labels = region_labels(sig_any_time);
    Significant_Regions.(contrast_names{c}).count = sum(sig_any_time);
    
    fprintf('\n%s: %d regions significantly different from zero (p < %.3f, Bonferroni corrected)\n', ...
        contrast_names{c}, sum(sig_any_time), Alpha);
end

%% 9. Create Summary Brain Map (Optional - requires CANLab tools)

% This creates a brain map showing region importance
% Uncomment if you want to create actual brain visualizations

 for c = 1:contrasts
     importance_vec = Beta_importance_overall(c, :)';
     
    % Create a new atlas object with importance values
     atlas_weighted = atlas_obj;
     atlas_weighted.probability_maps = importance_vec;
    
     % Visualize
     figure;
     montage(atlas_weighted, 'regioncenters');
     title(sprintf('Region Importance: %s', strrep(contrast_names{c}, '_', ' ')));
     
     saveas(gcf, sprintf('Brain_Map_%s.png', contrast_names{c}));
 end

%% 10. Save All Results

save('Region_Analysis_Results.mat', 'Top_Regions', 'Significant_Regions', ...
     'Beta_mean', 'Beta_std', 'Beta_se', 'Beta_consistency', ...
     'Beta_importance_overall', 'Beta_importance_peak', ...
     'Region_Consistency', 'multi_contrast_regions_sorted', ...
     'region_labels', 'contrast_names');

% Export top regions to CSV for each contrast
for c = 1:contrasts
    T = table();
    T.Rank = (1:top_n)';
    T.Region = Top_Regions.(contrast_names{c}).labels;
    T.Importance = Top_Regions.(contrast_names{c}).importance';
    
    writetable(T, sprintf('Top_Regions_%s.csv', contrast_names{c}));
end

fprintf('\n=== Analysis Complete ===\n');
fprintf('Results saved to Region_Analysis_Results.mat\n');
fprintf('Top regions per contrast exported to CSV files\n');
fprintf('Figures saved as PNG files\n');