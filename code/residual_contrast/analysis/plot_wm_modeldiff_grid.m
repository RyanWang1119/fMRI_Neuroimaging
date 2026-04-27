function plot_wm_modeldiff_grid(D, modelNames, fieldName, yLabel, saveFile, nPerm, TR)
% plot_wm_modeldiff_grid
% 6 WM contrasts x 3 model comparisons
%
% columns:
%   1 cHRF vs sHRF
%   2 cHRFderiv vs sHRF
%   3 cHRF vs cHRFderiv

if nargin < 6 || isempty(nPerm)
    nPerm = 5000;
end
if nargin < 7 || isempty(TR)
    TR = 0.72;
end

pairIdx = [1 3; 2 3; 1 2];
pairTitles = {'cHRF vs sHRF', 'cHRFderiv vs sHRF', 'cHRF vs cHRFderiv'};

contrast_labels = D{1}.contrast_labels;
K = numel(contrast_labels);

% assume [S x T x K]
T = size(D{1}.(fieldName), 2);
time = (0:T-1) * TR;

R = cell(K, 3);
allAbs = [];

for k = 1:K
    for c = 1:3
        a = pairIdx(c,1);
        b = pairIdx(c,2);

        X = squeeze(D{a}.(fieldName)(:,:,k));   % [S x T]
        Y = squeeze(D{b}.(fieldName)(:,:,k));   % [S x T]

        R{k,c} = cluster_perm_1d_paired(X, Y, nPerm, 0.05);
        allAbs = [allAbs, abs(R{k,c}.meanDiff)]; %#ok<AGROW>
    end
end

ymax = max(allAbs, [], 'omitnan');
if isempty(ymax) || ~isfinite(ymax) || ymax == 0
    ymax = 1;
end
ylim_use = 1.10 * [-ymax ymax];

fig = figure('Color','w','Position',[100 80 1050 1350]);
tiledlayout(K, 3, 'Padding','compact', 'TileSpacing','compact');

for k = 1:K
    for c = 1:3
        nexttile; hold on
        out = R{k,c};

        yl = ylim_use;
        yTopText = yl(2) - 0.12 * range(yl);
        yBotText = yl(1) + 0.12 * range(yl);

        % significant positive clusters
        for j = 1:numel(out.posClusters)
            cl = out.posClusters(j);
            if cl.p < 0.05
                patch([time(cl.start) time(cl.stop) time(cl.stop) time(cl.start)], ...
                      [yl(1) yl(1) yl(2) yl(2)], ...
                      [0.90 0.96 0.90], ...
                      'EdgeColor','none', 'FaceAlpha',0.8);
                text(mean(time(cl.start:cl.stop)), yTopText, p_text(cl.p), ...
                    'HorizontalAlignment','center', 'FontSize',7, 'Color',[0.15 0.35 0.15]);
            end
        end

        % significant negative clusters
        for j = 1:numel(out.negClusters)
            cl = out.negClusters(j);
            if cl.p < 0.05
                patch([time(cl.start) time(cl.stop) time(cl.stop) time(cl.start)], ...
                      [yl(1) yl(1) yl(2) yl(2)], ...
                      [0.97 0.90 0.90], ...
                      'EdgeColor','none', 'FaceAlpha',0.8);
                text(mean(time(cl.start:cl.stop)), yBotText, p_text(cl.p), ...
                    'HorizontalAlignment','center', 'FontSize',7, 'Color',[0.45 0.20 0.20]);
            end
        end

        plot(time, out.meanDiff, 'k-', 'LineWidth', 1.2);
        yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);

        xlim([time(1) time(end)]);
        ylim(ylim_use);
        grid off
        box off
        set(gca, 'FontSize', 8, 'LineWidth', 0.8);

        if k == 1
            title(pairTitles{c}, 'FontSize', 10, 'FontWeight', 'bold');
        end

        if c == 1
            ylabel({format_wm_contrast_label(contrast_labels{k}); yLabel}, 'FontSize', 8);
        end

        if k == K
            xlabel('Seconds', 'FontSize', 8);
        end
    end
end

exportgraphics(fig, saveFile, 'Resolution', 300);
end

function s = p_text(p)
if p <= 0.001
    s = 'p<0.001';
else
    s = sprintf('p=%.3f', p);
end
end

function s = format_wm_contrast_label(lbl)

lbl = char(lbl);

if contains(lbl, 'Body_LoadDiff_vs_Face_LoadDiff')
    s = 'Body vs Face';
elseif contains(lbl, 'Body_LoadDiff_vs_Place_LoadDiff')
    s = 'Body vs Place';
elseif contains(lbl, 'Body_LoadDiff_vs_Tool_LoadDiff')
    s = 'Body vs Tool';
elseif contains(lbl, 'Face_LoadDiff_vs_Place_LoadDiff')
    s = 'Face vs Place';
elseif contains(lbl, 'Face_LoadDiff_vs_Tool_LoadDiff')
    s = 'Face vs Tool';
elseif contains(lbl, 'Place_LoadDiff_vs_Tool_LoadDiff')
    s = 'Place vs Tool';
else
    s = strrep(lbl, '_LoadDiff', '');
    s = strrep(s, '_vs_', ' vs ');
    s = strrep(s, '_', ' ');
end

end