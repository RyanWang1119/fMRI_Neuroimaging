%% plot_WM_cluster_figures.m
% Assumes these are already in memory:
%   Tstab
%   results_all
%   stability_all
%
% Tstab columns expected:
% Model, Contrast, Cluster, Size, WithinConsensus, BetweenConsensus,
% ConsensusGap, MeanRegionStability, JaccardMean, JaccardSE, OverallARI

model_names    = {'sHRF','cHRFderiv','cHRF'};
contrast_names = {'body','faces','places','tools'};

%% Figure 1: cluster-size proportions across models
fig1 = figure('Color','w','Name','WM cluster size proportions');
tiledlayout(2,2, 'Padding','compact', 'TileSpacing','compact');

for c = 1:numel(contrast_names)
    nexttile;
    P = nan(numel(model_names), 2);   % [dominant minority]
    N = nan(numel(model_names), 2);   % raw sizes

    for m = 1:numel(model_names)
        sub = get_two_cluster_rows(Tstab, model_names{m}, contrast_names{c});
        [~, ord] = sort(sub.Size, 'descend');
        sub = sub(ord, :);  % dominant first, minority second

        N(m, :) = [sub.Size(1), sub.Size(2)];
        P(m, :) = N(m, :) ./ sum(N(m, :));
    end

    b = bar(P, 'stacked', 'LineWidth', 1);
    xticks(1:numel(model_names));
    xticklabels(model_names);
    ylim([0 1]);
    ylabel('Proportion of regions');
    title(sprintf('%s', contrast_names{c}));
    grid on;

    if c == 1
        legend({'Dominant cluster','Minority cluster'}, 'Location', 'southoutside');
    end

    % annotate raw sizes
    for m = 1:numel(model_names)
        text(m, 0.5, sprintf('%d / %d', N(m,1), N(m,2)), ...
            'HorizontalAlignment','center', 'FontSize',9, 'Color','k');
    end
end

sgtitle('WM region clustering: dominant vs minority cluster size proportions');


%% Figure 2: heatmaps of minority-cluster stability and overall ARI
minority_jaccard = nan(numel(contrast_names), numel(model_names));
minority_gap     = nan(numel(contrast_names), numel(model_names));
overall_ari      = nan(numel(contrast_names), numel(model_names));

for c = 1:numel(contrast_names)
    for m = 1:numel(model_names)
        sub = get_two_cluster_rows(Tstab, model_names{m}, contrast_names{c});
        [~, ord] = sort(sub.Size, 'descend');
        sub = sub(ord, :); % dominant first, minority second

        minority_jaccard(c,m) = sub.JaccardMean(2);
        minority_gap(c,m)     = sub.ConsensusGap(2);
        overall_ari(c,m)      = sub.OverallARI(1);  % same for both rows
    end
end

fig2 = figure('Color','w','Name','WM stability heatmaps');
tiledlayout(1,3, 'Padding','compact', 'TileSpacing','compact');

nexttile;
plot_heatmap(minority_jaccard, model_names, contrast_names, ...
    'Minority-cluster JaccardMean', [0 1]);

nexttile;
plot_heatmap(minority_gap, model_names, contrast_names, ...
    'Minority-cluster ConsensusGap', []);

nexttile;
plot_heatmap(overall_ari, model_names, contrast_names, ...
    'Overall ARI', [0 1]);

sgtitle('WM region clustering stability across HRF models');


%% Figure 3: cluster-average residual trajectories (4 x 3 panels)
fig3 = figure('Color','w','Name','WM cluster-average trajectories');
tiledlayout(numel(contrast_names), numel(model_names), ...
    'Padding','compact', 'TileSpacing','compact');

for c = 1:numel(contrast_names)
    for m = 1:numel(model_names)
        nexttile;

        res  = results_all.(model_names{m}){c};
        stab = stability_all.(model_names{m}){c};

        counts = res.cluster_counts(:);
        [~, ord] = sort(counts, 'descend');   % dominant first
        dom_k = ord(1);
        min_k = ord(2);

        plot(res.cluster_means(dom_k, :), 'LineWidth', 2); hold on;
        plot(res.cluster_means(min_k, :), 'LineWidth', 2);
        yline(0, '--', 'LineWidth', 0.8);
        grid on;

        if c == numel(contrast_names)
            xlabel('Time');
        end
        if m == 1
            ylabel(sprintf('%s\nMean residual', contrast_names{c}));
        end

        ct = stab.cluster_table;
        [~, ord_ct] = sort(ct.Size, 'descend');
        ct = ct(ord_ct, :);

        title(sprintf('%s\nsizes: %d/%d | minJ=%.2f | ARI=%.2f', ...
            model_names{m}, ct.Size(1), ct.Size(2), ...
            ct.JaccardMean(2), stab.overall_ari_mean), ...
            'FontSize', 10);

        if c == 1 && m == 1
            legend({'Dominant cluster','Minority cluster','Zero'}, ...
                'Location','best');
        end
    end
end

sgtitle('WM cluster-average residual trajectories by contrast and HRF model');


%% ---------- helper functions ----------

function sub = get_two_cluster_rows(Tstab, model_name, contrast_name)
    model_col    = string(Tstab.Model);
    contrast_col = string(Tstab.Contrast);

    idx = model_col == string(model_name) & contrast_col == string(contrast_name);
    sub = Tstab(idx, :);

    if height(sub) < 2
        error('Need at least 2 cluster rows for %s %s.', model_name, contrast_name);
    end

    % If somehow more than 2 rows are present, keep the two largest
    if height(sub) > 2
        [~, ord] = sort(sub.Size, 'descend');
        sub = sub(ord(1:2), :);
    end
end

function plot_heatmap(M, xlabels_in, ylabels_in, ttl, clim_in)
    imagesc(M);
    axis tight;
    xticks(1:numel(xlabels_in));
    xticklabels(xlabels_in);
    yticks(1:numel(ylabels_in));
    yticklabels(ylabels_in);
    xlabel('Model');
    ylabel('Contrast');
    title(ttl);
    colorbar;

    if ~isempty(clim_in)
        clim(clim_in);
    end

    for i = 1:size(M,1)
        for j = 1:size(M,2)
            if ~isnan(M(i,j))
                text(j, i, sprintf('%.2f', M(i,j)), ...
                    'HorizontalAlignment','center', ...
                    'FontSize', 9, 'Color', 'w', 'FontWeight', 'bold');
            end
        end
    end
end