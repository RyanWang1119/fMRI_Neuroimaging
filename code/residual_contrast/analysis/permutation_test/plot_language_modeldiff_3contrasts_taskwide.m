function plot_language_modeldiff_3contrasts_taskwide(D, fieldName, saveFile, nPerm, TR)

if nargin < 4 || isempty(nPerm), nPerm = 5000; end
if nargin < 5 || isempty(TR), TR = 0.72; end

time = (0:40) * TR;
contrastIdx = 1:3;

rowLabels = { ...
    'Present Phase: Story vs. Math', ...
    'Question Phase: Story vs. Math', ...
    'Response Phase: Story vs. Math'};

pairIdx = [1 3; 2 3; 1 2];
pairTitles = {'cHRF vs sHRF', 'cHRFderiv vs sHRF', 'cHRF vs cHRFderiv'};

R = cell(3,3);
allAbs = [];

for c = 1:3
    a = pairIdx(c,1);
    b = pairIdx(c,2);

    Xcell = cell(3,1);
    Ycell = cell(3,1);
    for r = 1:3
        k = contrastIdx(r);
        Xcell{r} = squeeze(D{a}.(fieldName)(:,:,k));
        Ycell{r} = squeeze(D{b}.(fieldName)(:,:,k));
    end

    tmp = cluster_perm_taskwide_paired(Xcell, Ycell, nPerm, 0.05);

    for r = 1:3
        R{r,c} = tmp(r);
        allAbs = [allAbs, abs(tmp(r).meanDiff)]; %#ok<AGROW>
    end
end

ymax = max(allAbs, [], 'omitnan');
if isempty(ymax) || ~isfinite(ymax) || ymax == 0
    ymax = 1;
end
ylim_use = 1.15 * [-ymax ymax];

fig = figure('Color','w','Position',[100 100 1050 780]);
tiledlayout(3,3,'Padding','compact','TileSpacing','compact');

for r = 1:3
    for c = 1:3
        nexttile; hold on
        out = R{r,c};

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
                    'HorizontalAlignment','center', 'FontSize',7, 'Color',[0.2 0.35 0.2]);
            end
        end

        for j = 1:numel(out.negClusters)
            cl = out.negClusters(j);
            if cl.p < 0.05
                patch([time(cl.start) time(cl.stop) time(cl.stop) time(cl.start)], ...
                      [yl(1) yl(1) yl(2) yl(2)], ...
                      [0.97 0.90 0.90], 'EdgeColor','none', 'FaceAlpha',0.8);
                text(mean(time(cl.start:cl.stop)), yBotText, p_text(cl.p), ...
                    'HorizontalAlignment','center', 'FontSize',7, 'Color',[0.45 0.2 0.2]);
            end
        end

        plot(time, out.meanDiff, 'k-', 'LineWidth', 1.3);
        yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);

        xlim([time(1) time(end)]);
        ylim(ylim_use);
        grid off
        box off
        set(gca, 'FontSize', 8, 'LineWidth', 0.8);

        if r == 1
            title(pairTitles{c}, 'FontSize', 13, 'FontWeight', 'bold');
        end
        if c == 1
            ylabel(rowLabels{r}, 'FontSize', 12);
        end
        if r == 3
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