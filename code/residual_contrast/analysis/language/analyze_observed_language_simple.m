function summaryTbl = analyze_observed_language_simple(baseDir, outdir)
% analyze_observed_language_simple
% Local downstream analysis for:
%   observed_language_simple_chrf.mat
%   observed_language_simple_chrfderiv.mat
%   observed_language_simple_shrf.mat

if nargin < 1 || isempty(baseDir)
    baseDir = pwd;
end
if nargin < 2 || isempty(outdir)
    outdir = fullfile(baseDir, 'language_simple_analysis');
end
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

obsFiles = struct( ...
    'cHRF',      fullfile(baseDir, 'observed_language_simple_chrf.mat'), ...
    'cHRFderiv', fullfile(baseDir, 'observed_language_simple_chrfderiv.mat'), ...
    'sHRF',      fullfile(baseDir, 'observed_language_simple_shrf.mat'));

modelNames = fieldnames(obsFiles);
nModels = numel(modelNames);

data = cell(nModels,1);
for m = 1:nModels
    f = obsFiles.(modelNames{m});
    if ~isfile(f)
        error('Missing file: %s', f);
    end
    S = load(f);
    data{m} = check_and_standardize(S, modelNames{m});
end

% check contrast labels match across models
refLabels = data{1}.contrast_labels(:);
for m = 2:nModels
    if ~isequal(refLabels, data{m}.contrast_labels(:))
        error('contrast_labels do not match across models.');
    end
end

contrast_labels = refLabels;
T = size(data{1}.TA_mean, 1);
time_s = (0:T-1) * 0.72;

summaryTbl = build_summary_table(data, modelNames, contrast_labels, time_s);
writetable(summaryTbl, fullfile(outdir, 'language_simple_summary.csv'));

plot_metric_grid(data, modelNames, contrast_labels, time_s, ...
    'TA_mean', 'Temporal Accuracy', true, fullfile(outdir, 'language_simple_TA.png'));

hasPMZ = all(cellfun(@(x) isfield(x,'PMZ_mean') && ~isempty(x.PMZ_mean), data));
if hasPMZ
    plot_metric_grid(data, modelNames, contrast_labels, time_s, ...
        'PMZ_mean', 'PMZ', false, fullfile(outdir, 'language_simple_PMZ.png'));
end

make_haufe_heatmaps(data, modelNames, contrast_labels, time_s, outdir);

save(fullfile(outdir, 'language_simple_summary.mat'), ...
    'summaryTbl', 'contrast_labels', 'time_s', 'modelNames');

fprintf('Saved local LANGUAGE analysis to:\n%s\n', outdir);

end

% =========================
% helpers
% =========================
function D = check_and_standardize(S, modelName)

req = {'TA_mean','Haufe_mean','contrast_labels'};
for i = 1:numel(req)
    if ~isfield(S, req{i})
        error('%s missing field: %s', modelName, req{i});
    end
end

D = struct();
D.TA_mean = S.TA_mean;
D.Haufe_mean = S.Haufe_mean;
D.contrast_labels = cellstr(S.contrast_labels);
if isfield(S, 'PMZ_mean')
    D.PMZ_mean = S.PMZ_mean;
else
    D.PMZ_mean = [];
end

% expected shapes
if ndims(D.TA_mean) ~= 2
    error('%s: TA_mean must be 2D.', modelName);
end
if size(D.TA_mean,1) ~= 41 && size(D.TA_mean,2) == 41
    D.TA_mean = D.TA_mean.';
end
if size(D.TA_mean,1) ~= 41
    error('%s: TA_mean should have 41 time points.', modelName);
end

if ~isempty(D.PMZ_mean)
    if size(D.PMZ_mean,1) ~= 41 && size(D.PMZ_mean,2) == 41
        D.PMZ_mean = D.PMZ_mean.';
    end
    if ~isequal(size(D.PMZ_mean), size(D.TA_mean))
        error('%s: PMZ_mean size mismatch.', modelName);
    end
end

if ndims(D.Haufe_mean) ~= 3
    error('%s: Haufe_mean must be 3D.', modelName);
end

sz = size(D.Haufe_mean);
K = size(D.TA_mean,2);

if isequal(sz, [489 41 K])
    % ok
elseif isequal(sz, [41 489 K])
    D.Haufe_mean = permute(D.Haufe_mean, [2 1 3]);
elseif isequal(sz, [489 K 41])
    D.Haufe_mean = permute(D.Haufe_mean, [1 3 2]);
else
    error('%s: unrecognized Haufe_mean size: [%s]', modelName, num2str(sz));
end

if numel(D.contrast_labels) ~= K
    error('%s: number of contrast_labels does not match TA_mean.', modelName);
end

end

function summaryTbl = build_summary_table(data, modelNames, contrast_labels, time_s)

rows = [];
for m = 1:numel(data)
    D = data{m};
    for k = 1:numel(contrast_labels)
        [peakTA, idxTA] = max(D.TA_mean(:,k));
        meanTA = mean(D.TA_mean(:,k), 'omitnan');

        if isfield(D,'PMZ_mean') && ~isempty(D.PMZ_mean)
            [peakPMZ, idxPMZ] = max(D.PMZ_mean(:,k));
            meanPMZ = mean(D.PMZ_mean(:,k), 'omitnan');
        else
            peakPMZ = NaN;
            idxPMZ = NaN;
            meanPMZ = NaN;
        end

        row = table( ...
            string(modelNames{m}), ...
            string(contrast_labels{k}), ...
            peakTA, idxTA, time_s(idxTA), meanTA, ...
            peakPMZ, idxPMZ, iff(~isnan(idxPMZ), time_s(idxPMZ), NaN), meanPMZ, ...
            'VariableNames', {'Model','Contrast','PeakTA','PeakTA_idx','PeakTA_sec','MeanTA', ...
                              'PeakPMZ','PeakPMZ_idx','PeakPMZ_sec','MeanPMZ'});
        rows = [rows; row];
    end
end

summaryTbl = rows;
end

function plot_metric_grid(data, modelNames, contrast_labels, time_s, fieldName, ylab, addChance, saveFile)

K = numel(contrast_labels);
nRow = 3;
nCol = 3;

fig = figure('Color','w','Position',[100 100 1500 1000]);
tiledlayout(nRow, nCol, 'Padding','compact', 'TileSpacing','compact');

for k = 1:K
    nexttile;
    hold on
    for m = 1:numel(data)
        plot(time_s, data{m}.(fieldName)(:,k), 'LineWidth', 1.8);
    end
    if addChance
        yline(0.5, 'k--', 'LineWidth', 1);
    end
    title(strrep(contrast_labels{k}, '_', '\_'), 'Interpreter','tex', 'FontSize',10);
    xlabel('Time (s)');
    ylabel(ylab);
    xlim([time_s(1) time_s(end)]);
    grid on
end

legend(modelNames, 'Location','southoutside', 'Orientation','horizontal');
exportgraphics(fig, saveFile, 'Resolution', 300);
close(fig);
end

function make_haufe_heatmaps(data, modelNames, contrast_labels, time_s, outdir)

allVals = [];
for m = 1:numel(data)
    allVals = [allVals; data{m}.Haufe_mean(:)];
end
absLim = prctile(abs(allVals), 99);
if absLim <= 0
    absLim = max(abs(allVals));
end
cax = [-absLim absLim];

for k = 1:numel(contrast_labels)
    fig = figure('Color','w','Position',[100 100 1500 900]);
    tiledlayout(numel(data), 1, 'Padding','compact', 'TileSpacing','compact');

    for m = 1:numel(data)
        nexttile;
        imagesc(time_s, 1:size(data{m}.Haufe_mean,1), data{m}.Haufe_mean(:,:,k));
        set(gca, 'YDir', 'normal');
        clim(cax);
        colorbar;
        xlabel('Time (s)');
        ylabel('Parcel');
        title(sprintf('%s | %s', modelNames{m}, contrast_labels{k}), ...
            'Interpreter','none');
    end

    exportgraphics(fig, fullfile(outdir, sprintf('haufe_%02d_%s.png', ...
        k, sanitize_filename(contrast_labels{k}))), 'Resolution', 300);
    close(fig);
end
end

function out = sanitize_filename(s)
out = regexprep(s, '[^a-zA-Z0-9_-]', '_');
end

function y = iff(cond, a, b)
if cond
    y = a;
else
    y = b;
end
end