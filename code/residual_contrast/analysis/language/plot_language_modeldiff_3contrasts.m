function plot_language_modeldiff_3contrasts(D, modelNames, fieldName, yLabel, figTitle, saveFile, nPerm)

time = (0:40)*0.72;
contrastIdx = 1:3;

rowLabels = { ...
    'Present Phase: Story vs. Math', ...
    'Question Phase: Story vs. Math', ...
    'Response Phase: Story vs. Math'};

pairIdx = [1 3; 2 3; 1 2];
pairTitles = { ...
    'cHRF vs sHRF', ...
    'cHRFderiv vs sHRF', ...
    'cHRF vs cHRFderiv'};

% precompute all results first
R = cell(3,3);
allAbs = [];

for r = 1:3
    k = contrastIdx(r);
    for c = 1:3
        a = pairIdx(c,1);
        b = pairIdx(c,2);

        X = squeeze(D{a}.(fieldName)(:,:,k));   % [S x T]
        Y = squeeze(D{b}.(fieldName)(:,:,k));   % [S x T]

        R{r,c} = cluster_perm_1d_paired(X, Y, nPerm, 0.05);
        allAbs = [allAbs, abs(R{r,c}.meanDiff)];
    end
end

ymax = max(allAbs);
if isempty(ymax) || ymax == 0
    ymax = 1;
end
ylim_use = 1.15 * [-ymax ymax];

fig = figure('Color','w','Position',[100 100 1050 780]);
tiledlayout(3,3,'Padding','compact','TileSpacing','compact');

for r = 1:3
    for c = 1:3
        nexttile; hold on
        out = R{r,c};

        % background shading for significant clusters
        yl = ylim_use;

        for j = 1:numel(out.posClusters)
            cl = out.posClusters(j);
            if cl.p < 0.05
                patch([time(cl.start) time(cl.stop) time(cl.stop) time(cl.start)], ...
                      [yl(1) yl(1) yl(2) yl(2)], ...
                      [0.90 0.96 0.90], 'EdgeColor','none', 'FaceAlpha',0.8);
                text(mean(time(cl.start:cl.stop)), yl(2)*0.82, p_text(cl.p), ...
                    'HorizontalAlignment','center', 'FontSize',7, 'Color',[0.2 0.35 0.2]);
            end
        end

        for j = 1:numel(out.negClusters)
            cl = out.negClusters(j);
            if cl.p < 0.05
                patch([time(cl.start) time(cl.stop) time(cl.stop) time(cl.start)], ...
                      [yl(1) yl(1) yl(2) yl(2)], ...
                      [0.97 0.90 0.90], 'EdgeColor','none', 'FaceAlpha',0.8);
                text(mean(time(cl.start:cl.stop)), yl(1)*0.82, p_text(cl.p), ...
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
            title(pairTitles{c}, 'FontSize', 10, 'FontWeight', 'bold');
        end
        if c == 1
            ylabel(rowLabels{r}, 'FontSize', 8);
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