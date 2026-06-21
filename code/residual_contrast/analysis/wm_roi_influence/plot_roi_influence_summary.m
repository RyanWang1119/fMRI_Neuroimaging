function plot_roi_influence_summary(summaryTable, outerResults, cfg, outDir)
%PLOT_ROI_INFLUENCE_SUMMARY Create required ROI influence summary plots.

if nargin < 4 || isempty(outDir)
    outDir = pwd;
end
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

TAFC = arrayfun(@(x) x.TAFC, outerResults(:));

f = figure('Color', 'w', 'Visible', 'off');
histogram(TAFC, 'BinWidth', 0.025);
xlabel('Held-out subject-level TAFC');
ylabel('Outer splits');
title('Outer-split TAFC distribution');
xline(0.5, 'k--', 'Chance');
save_figure(f, fullfile(outDir, 'outer_split_tafc_distribution.png'));

plot_ranked_bar(summaryTable, 'mean_model_reliance_tafc', ...
    'Mean held-out model reliance (delta TAFC)', ...
    'parcel_reliance_ranking.png', cfg, outDir);

plot_ranked_bar(summaryTable, 'selection_frequency', ...
    'Elastic-net selection frequency', ...
    'selection_frequency_ranking.png', cfg, outDir);

f = figure('Color', 'w', 'Visible', 'off');
scatter(summaryTable.selection_frequency, summaryTable.mean_model_reliance_tafc, ...
    28, double(summaryTable.high_confidence_influential_candidate), 'filled');
xlabel('Elastic-net selection frequency');
ylabel('Mean held-out model reliance (delta TAFC)');
title('Selection frequency versus held-out model reliance');
xline(cfg.thresholds.minSelectionFrequency, 'k--');
yline(cfg.thresholds.minMeanReliance, 'k--');
grid on;
save_figure(f, fullfile(outDir, 'selection_vs_model_reliance.png'));
end

function plot_ranked_bar(T, metricName, ylab, fileName, cfg, outDir)
vals = T.(metricName);
[~, ord] = sort(vals, 'descend', 'MissingPlacement', 'last');
topN = min(cfg.plots.topN, height(T));
ord = ord(1:topN);

labels = T.roi_label(ord);
if isstring(labels) || iscategorical(labels)
    labels = cellstr(labels);
end
for i = 1:numel(labels)
    labels{i} = sprintf('%d: %s', T.roi_index(ord(i)), labels{i});
end

f = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 1200 600]);
bar(vals(ord));
ylabel(ylab);
title(sprintf('Top %d parcels by %s', topN, strrep(metricName, '_', ' ')));
set(gca, 'XTick', 1:topN, 'XTickLabel', labels, 'XTickLabelRotation', 60);
grid on;
save_figure(f, fullfile(outDir, fileName));
end

function save_figure(f, fileName)
try
    exportgraphics(f, fileName, 'Resolution', 300);
catch
    saveas(f, fileName);
end
close(f);
end
