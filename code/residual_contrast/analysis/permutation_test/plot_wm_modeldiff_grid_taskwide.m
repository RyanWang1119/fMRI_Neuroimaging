function plot_wm_modeldiff_grid_taskwide(D, fieldName, saveFile, nPerm, TR)

if nargin < 4 || isempty(nPerm), nPerm = 5000; end
if nargin < 5 || isempty(TR), TR = 0.72; end

pairIdx = [1 3; 2 3; 1 2];
pairTitles = {'cHRF vs sHRF', 'cHRFderiv vs sHRF', 'cHRF vs cHRFderiv'};

contrast_labels = D{1}.contrast_labels;
K = numel(contrast_labels);
T = size(D{1}.(fieldName), 2);
time = (0:T-1) * TR;

R = cell(K,3);
allAbs = [];

for c = 1:3
    a = pairIdx(c,1);
    b = pairIdx(c,2);

    Xcell = cell(K,1);
    Ycell = cell(K,1);
    for k = 1:K
        Xcell{k} = squeeze(D{a}.(fieldName)(:,:,k));
        Ycell{k} = squeeze(D{b}.(fieldName)(:,:,k));
    end

    tmp = cluster_perm_taskwide_paired(Xcell, Ycell, nPerm, 0.05);

    for k = 1:K
        R{k,c} = tmp(k);
        allAbs = [allAbs, abs(tmp(k).meanDiff)]; %#ok<AGROW>
    end
end

ymax = max(allAbs, [], 'omitnan');
if isempty(ymax) || ~isfinite(ymax) || ymax == 0
    ymax = 1;
end
ylim_use = 1.10 * [-ymax ymax];

fig = figure('Color','w','Position',[100 80 1050 1350]);
tiledlayout(K,3,'Padding','compact','TileSpacing','compact');

for k = 1:K
    for c = 1:3
        nexttile; hold on
        out = R{k,c};

        yl = ylim_use;
        yTopText = yl(2) - 0.12 * range(yl);
        yBotText = yl(1) + 0.12 * range(yl);

        for j = 1:numel(out.posClusters)
            cl = out.posClusters(j);
            if cl.p < 0.05
                patch([time(cl.start) time(cl.stop) time(cl.stop) time(cl.start)], ...
                      [yl(1) yl(1) yl(2) yl(2)], ...
                      [0.90 0.96 0.90], 'EdgeColor','none', 'FaceAlpha',0.8);
                text(mean(time(cl.start:cl.stop)), yTopText, p_text(cl.p), ...
                    'HorizontalAlignment','center', 'FontSize',7, 'Color',[0.15 0.35 0.15]);
            end
        end

        for j = 1:numel(out.negClusters)
            cl = out.negClusters(j);
            if cl.p < 0.05
                patch([time(cl.start) time(cl.stop) time(cl.stop) time(cl.start)], ...
                      [yl(1) yl(1) yl(2) yl(2)], ...
                      [0.97 0.90 0.90], 'EdgeColor','none', 'FaceAlpha',0.8);
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
            title(pairTitles{c}, 'FontSize', 13, 'FontWeight', 'bold');
        end
        if c == 1
            ylabel(format_wm_contrast_label(contrast_labels{k}), 'FontSize', 12);
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
lbl = strrep(lbl, 'Faces',  'Face');
lbl = strrep(lbl, 'Places', 'Place');
lbl = strrep(lbl, 'Tools',  'Tool');
lbl = strrep(lbl, '_LoadDiff', '');
lbl = strrep(lbl, '_vs_', '|');
lbl = strrep(lbl, '_', '');

switch lower(lbl)
    case 'body|face'
        s = 'Body vs Face';
    case 'body|place'
        s = 'Body vs Place';
    case 'body|tool'
        s = 'Body vs Tool';
    case 'face|place'
        s = 'Face vs Place';
    case 'face|tool'
        s = 'Face vs Tool';
    case 'place|tool'
        s = 'Place vs Tool';
    otherwise
        s = strrep(lbl, '|', ' vs ');
end
end