function summaryTable = summarize_roi_influence(outerResults, signFlipResults, roiLabels, cfg)
%SUMMARIZE_ROI_INFLUENCE Build one-row-per-parcel candidate summary.
%
% Outputs are intentionally separated:
%   1. Haufe pattern: task-linked residual signal map.
%   2. Elastic-net selection frequency: stability of sparse multivariable selection.
%   3. Held-out model reliance: dependence of this fitted decoder on that parcel.
%
% These thresholds define a reproducibility-based candidate set; they are
% not causal, mechanistic, exclusive, or definitive biological region labels.

nRepeats = numel(outerResults);
nROI = numel(outerResults(1).betaStandardized);

Hraw = nan(nROI, nRepeats);
Bstd = nan(nROI, nRepeats);
Sel = false(nROI, nRepeats);

nSwaps = size(outerResults(1).modelRelianceTAFC, 2);
Reliance = zeros(nROI, nSwaps, nRepeats, 'single');

for r = 1:nRepeats
    Hraw(:, r) = outerResults(r).haufePatternRaw(:);
    Bstd(:, r) = outerResults(r).betaStandardized(:);
    Sel(:, r) = outerResults(r).selectionIndicator(:);
    Reliance(:, :, r) = outerResults(r).modelRelianceTAFC;
end

selection_frequency = mean(Sel, 2);
mean_haufe_raw = mean(Hraw, 2, 'omitnan');
median_haufe_raw = median(Hraw, 2, 'omitnan');
haufe_sign_consistency = sign_consistency(Hraw);
beta_sign_consistency = sign_consistency_selected(Bstd, Sel);

R2 = reshape(Reliance, nROI, []);
mean_model_reliance_tafc = mean(double(R2), 2, 'omitnan');
median_model_reliance_tafc = median(double(R2), 2, 'omitnan');

stable_predictive_candidate = ...
    selection_frequency >= cfg.thresholds.minSelectionFrequency & ...
    beta_sign_consistency >= cfg.thresholds.minSignConsistency & ...
    mean_model_reliance_tafc > cfg.thresholds.minMeanReliance;

population_supported = signFlipResults.pFWER < cfg.permutation.alpha;
if cfg.thresholds.requireFWER
    high_confidence_influential_candidate = stable_predictive_candidate & population_supported;
else
    high_confidence_influential_candidate = stable_predictive_candidate;
end

roiTable = normalize_roi_table(roiLabels, nROI);

summaryTable = table();
summaryTable.roi_index = roiTable.roi_index;
summaryTable.roi_label = roiTable.roi_label;
summaryTable.roi_network = roiTable.roi_network;
summaryTable.n_valid_subjects = signFlipResults.nValid(:);
summaryTable.mean_raw_difference = signFlipResults.meanRawDifference(:);
summaryTable.cohen_dz = signFlipResults.effectSizeCohenDz(:);
summaryTable.t_observed = signFlipResults.tObs(:);
summaryTable.p_fwer = signFlipResults.pFWER(:);
summaryTable.significant_fwer = signFlipResults.significantFWER(:);
summaryTable.mean_haufe_raw = mean_haufe_raw(:);
summaryTable.median_haufe_raw = median_haufe_raw(:);
summaryTable.haufe_sign_consistency = haufe_sign_consistency(:);
summaryTable.selection_frequency = selection_frequency(:);
summaryTable.beta_sign_consistency = beta_sign_consistency(:);
summaryTable.mean_model_reliance_tafc = mean_model_reliance_tafc(:);
summaryTable.median_model_reliance_tafc = median_model_reliance_tafc(:);
summaryTable.stable_predictive_candidate = stable_predictive_candidate(:);
summaryTable.population_supported = population_supported(:);
summaryTable.high_confidence_influential_candidate = high_confidence_influential_candidate(:);
end

function sc = sign_consistency(X)
n = size(X, 1);
sc = zeros(n, 1);
for i = 1:n
    x = X(i, :);
    x = x(isfinite(x) & x ~= 0);
    if isempty(x)
        sc(i) = 0;
    else
        sc(i) = max(mean(x > 0), mean(x < 0));
    end
end
end

function sc = sign_consistency_selected(B, Sel)
n = size(B, 1);
sc = zeros(n, 1);
for i = 1:n
    x = B(i, Sel(i, :) & isfinite(B(i, :)) & B(i, :) ~= 0);
    if isempty(x)
        sc(i) = 0;
    else
        sc(i) = max(mean(x > 0), mean(x < 0));
    end
end
end

function T = normalize_roi_table(roiLabels, nROI)
if istable(roiLabels)
    T = roiLabels;
else
    if isstring(roiLabels) || iscategorical(roiLabels)
        labels = cellstr(roiLabels(:));
    elseif iscell(roiLabels)
        labels = roiLabels(:);
    else
        labels = arrayfun(@(r) sprintf('Parcel_%03d', r), (1:nROI)', 'UniformOutput', false);
    end
    T = table((1:nROI)', labels(:), repmat({''}, nROI, 1), ...
        'VariableNames', {'roi_index', 'roi_label', 'roi_network'});
end

if height(T) ~= nROI
    error('roiLabels has %d rows; expected %d.', height(T), nROI);
end
if ~ismember('roi_index', T.Properties.VariableNames)
    T.roi_index = (1:nROI)';
end
if ~ismember('roi_label', T.Properties.VariableNames)
    T.roi_label = arrayfun(@(r) sprintf('Parcel_%03d', r), (1:nROI)', 'UniformOutput', false);
end
if ~ismember('roi_network', T.Properties.VariableNames)
    T.roi_network = repmat({''}, nROI, 1);
end
T = T(:, {'roi_index', 'roi_label', 'roi_network'});
end
