%% run_WM_region_overlap_analysis.m
% Requires:
%   results_all
%   stability_all
%   atlas_obj
%
% Uses either:
%   - raw minority cluster
%   - stability-aware substantive minority cluster

clearvars -except results_all stability_all atlas_obj

model_names = {'sHRF', 'cHRFderiv', 'cHRF'};
contrast_labels = {'body', 'faces', 'places', 'tools'};
nRegions = numel(atlas_obj.labels);

%% Settings
use_stability_aware = false;   % recommended
jaccard_thresh = 0.25;        % minimum cluster reproducibility
min_cluster_size = 5;         % minimum size for substantive minority

%% Run overlap analysis
overlap_results = struct();

for c = 1:numel(contrast_labels)
    contrast_name = contrast_labels{c};

    fprintf('\n=============================\n');
    fprintf('Contrast: %s\n', contrast_name);
    fprintf('=============================\n');

    membership = false(nRegions, numel(model_names));
    chosen_info = cell(numel(model_names), 1);

    for m = 1:numel(model_names)
        model_name = model_names{m};
        res  = results_all.(model_name){c};
        stab = stability_all.(model_name){c};

        if use_stability_aware
            [idx_orig, info] = get_substantive_minority_cluster(res, stab, ...
                'MinClusterSize', min_cluster_size, ...
                'JaccardThreshold', jaccard_thresh);
        else
            [idx_orig, info] = get_raw_minority_cluster(res);
        end

        membership(idx_orig, m) = true;
        chosen_info{m} = info;

        fprintf('%s: selected cluster %d, size=%d, Jaccard=%.3f\n', ...
            model_name, info.Cluster, info.Size, info.JaccardMean);
    end

    % Pairwise overlaps
    pairwise_table = compute_pairwise_overlap_table(membership, model_names);

    % Triple overlap
    triple_idx = find(all(membership, 2));
    triple_names = atlas_obj.labels(triple_idx);

    % Exact pattern table
    pattern_table = compute_pattern_table(membership, atlas_obj.labels, model_names);

    % Regions appearing in 2 or more models
    recur_idx = find(sum(membership, 2) >= 2);
    recur_names = atlas_obj.labels(recur_idx);
    recur_count = sum(membership(recur_idx, :), 2);

    recurrent_table = table(recur_idx(:), recur_names(:), recur_count(:), ...
        membership(recur_idx,1), membership(recur_idx,2), membership(recur_idx,3), ...
        'VariableNames', {'ParcelIndex','ParcelName','NModels','sHRF','cHRFderiv','cHRF'});

    overlap_results.(contrast_name).membership = membership;
    overlap_results.(contrast_name).chosen_info = chosen_info;
    overlap_results.(contrast_name).pairwise_table = pairwise_table;
    overlap_results.(contrast_name).triple_idx = triple_idx;
    overlap_results.(contrast_name).triple_names = triple_names;
    overlap_results.(contrast_name).pattern_table = pattern_table;
    overlap_results.(contrast_name).recurrent_table = recurrent_table;

    disp(pairwise_table)
    fprintf('Shared by all 3 models: %d regions\n', numel(triple_idx));
    if ~isempty(triple_idx)
        disp(table(triple_idx(:), triple_names(:), ...
            'VariableNames', {'ParcelIndex','ParcelName'}))
    end

    fprintf('Regions appearing in >=2 models: %d\n', height(recurrent_table));
end

%% Build summary table across contrasts
T_overlap_summary = summarize_overlap_results(overlap_results, contrast_labels);
disp(T_overlap_summary)

%% Plots
fig1 = plot_overlap_pattern_bars(overlap_results, contrast_labels);
fig2 = plot_pairwise_jaccard_heatmaps(overlap_results, contrast_labels);

%% Save
save('WM_region_overlap_results.mat', 'overlap_results', 'T_overlap_summary', '-v7.3');
writetable(T_overlap_summary, 'WM_region_overlap_summary.csv');

fprintf('\nSaved:\n');
fprintf('  WM_region_overlap_results.mat\n');
fprintf('  WM_region_overlap_summary.csv\n');

%% =========================
%% Helper functions
%% =========================

function [idx_orig, info] = get_raw_minority_cluster(res)
    counts = accumarray(res.cluster_labels, 1);
    [~, minority_cluster] = min(counts);

    local_idx = find(res.cluster_labels == minority_cluster);
    idx_orig = res.keep_idx(local_idx);

    info = table2struct(table(minority_cluster, numel(idx_orig), NaN, ...
        'VariableNames', {'Cluster','Size','JaccardMean'}));
end

function [idx_orig, info] = get_substantive_minority_cluster(res, stab, varargin)
    p = inputParser;
    addParameter(p, 'MinClusterSize', 5, @isscalar);
    addParameter(p, 'JaccardThreshold', 0.25, @isscalar);
    parse(p, varargin{:});

    min_size = p.Results.MinClusterSize;
    jac_thr = p.Results.JaccardThreshold;

    ct = stab.cluster_table;

    % sort by size descending
    [~, ord] = sort(ct.Size, 'descend');
    ct = ct(ord, :);

    % dominant cluster = largest
    dominant_cluster = ct.Cluster(1);

    % candidates = non-dominant clusters
    cand = ct(ct.Cluster ~= dominant_cluster, :);

    % keep substantive ones
    keep = cand.Size >= min_size & cand.JaccardMean >= jac_thr;

    if any(keep)
        cand2 = cand(keep, :);
        % choose the largest substantive non-dominant cluster
        [~, best] = max(cand2.Size);
        chosen_cluster = cand2.Cluster(best);
        chosen_size = cand2.Size(best);
        chosen_jac = cand2.JaccardMean(best);
    else
        % if no substantive cluster exists, return empty
        idx_orig = [];
        info = table2struct(table(NaN, 0, NaN, ...
            'VariableNames', {'Cluster','Size','JaccardMean'}));
        return
    end

    local_idx = find(res.cluster_labels == chosen_cluster);
    idx_orig = res.keep_idx(local_idx);

    info = table2struct(table(chosen_cluster, chosen_size, chosen_jac, ...
        'VariableNames', {'Cluster','Size','JaccardMean'}));
end

function Tpair = compute_pairwise_overlap_table(membership, model_names)
    pairs = nchoosek(1:numel(model_names), 2);
    rows = cell(size(pairs,1), 6);

    for i = 1:size(pairs,1)
        a = pairs(i,1);
        b = pairs(i,2);

        A = membership(:,a);
        B = membership(:,b);

        inter = sum(A & B);
        union = sum(A | B);

        if union == 0
            jac = NaN;
        else
            jac = inter / union;
        end

        rows(i,:) = {model_names{a}, model_names{b}, sum(A), sum(B), inter, jac};
    end

    Tpair = cell2table(rows, ...
        'VariableNames', {'Model1','Model2','N1','N2','Intersection','Jaccard'});
end

function Tpat = compute_pattern_table(membership, labels, model_names)

    % Build exact 3-bit pattern strings: 000, 100, 010, ..., 111
    pat_str = strings(size(membership,1), 1);
    for i = 1:size(membership,1)
        pat_str(i) = sprintf('%d%d%d', ...
            membership(i,1), membership(i,2), membership(i,3));
    end

    unique_patterns = ["000","100","010","001","110","101","011","111"];
    desc = { ...
        'none', ...
        [model_names{1} ' only'], ...
        [model_names{2} ' only'], ...
        [model_names{3} ' only'], ...
        [model_names{1} '+' model_names{2}], ...
        [model_names{1} '+' model_names{3}], ...
        [model_names{2} '+' model_names{3}], ...
        'all 3'};

    counts = zeros(numel(unique_patterns),1);
    examples = repmat({''}, numel(unique_patterns), 1);

    for i = 1:numel(unique_patterns)
        idx = find(pat_str == unique_patterns(i));
        counts(i) = numel(idx);

        if ~isempty(idx)
            nshow = min(5, numel(idx));
            examples{i} = strjoin(string(labels(idx(1:nshow))), '; ');
        end
    end

    Tpat = table(unique_patterns(:), desc(:), counts, examples, ...
        'VariableNames', {'Pattern','Description','Count','ExampleParcels'});
end

function T = summarize_overlap_results(overlap_results, contrast_labels)
    rows = {};

    for c = 1:numel(contrast_labels)
        nm = contrast_labels{c};
        R = overlap_results.(nm);

        membership = R.membership;

        n_s = sum(membership(:,1));
        n_d = sum(membership(:,2));
        n_c = sum(membership(:,3));

        inter_sd = sum(membership(:,1) & membership(:,2));
        inter_sc = sum(membership(:,1) & membership(:,3));
        inter_dc = sum(membership(:,2) & membership(:,3));
        inter_all = sum(all(membership,2));

        jac_sd = safe_jaccard(membership(:,1), membership(:,2));
        jac_sc = safe_jaccard(membership(:,1), membership(:,3));
        jac_dc = safe_jaccard(membership(:,2), membership(:,3));

        rows(end+1,:) = { ...
            nm, ...
            n_s, n_d, n_c, ...
            inter_sd, jac_sd, ...
            inter_sc, jac_sc, ...
            inter_dc, jac_dc, ...
            inter_all ...
            };
    end

    T = cell2table(rows, ...
        'VariableNames', { ...
        'Contrast', ...
        'NsHRF', 'NcHRFderiv', 'NcHRF', ...
        'Inter_sHRF_cHRFderiv', 'Jaccard_sHRF_cHRFderiv', ...
        'Inter_sHRF_cHRF', 'Jaccard_sHRF_cHRF', ...
        'Inter_cHRFderiv_cHRF', 'Jaccard_cHRFderiv_cHRF', ...
        'IntersectionAll3'});
end

function j = safe_jaccard(A, B)
    inter = sum(A & B);
    union = sum(A | B);
    if union == 0
        j = NaN;
    else
        j = inter / union;
    end
end

function fig = plot_overlap_pattern_bars(overlap_results, contrast_labels)
    fig = figure('Color','w', 'Name', 'WM overlap patterns');
    tiledlayout(2,2, 'Padding','compact', 'TileSpacing','compact');

    for c = 1:numel(contrast_labels)
        nexttile;
        Tpat = overlap_results.(contrast_labels{c}).pattern_table;

        % exclude 000 to focus on minority membership
        Tplot = Tpat(Tpat.Pattern ~= "000", :);

        bar(Tplot.Count);
        xticks(1:height(Tplot));
        xticklabels(Tplot.Description);
        xtickangle(45);
        ylabel('Number of regions');
        title(contrast_labels{c});
        grid on;
    end

    sgtitle('Exact minority-membership patterns across models');
end

function fig = plot_pairwise_jaccard_heatmaps(overlap_results, contrast_labels)
    fig = figure('Color','w', 'Name', 'WM pairwise overlap heatmaps');
    tiledlayout(1,3, 'Padding','compact', 'TileSpacing','compact');

    mats = zeros(numel(contrast_labels), 3);
    for c = 1:numel(contrast_labels)
        Tpair = overlap_results.(contrast_labels{c}).pairwise_table;
        mats(c,1) = Tpair.Jaccard(1); % sHRF vs cHRFderiv
        mats(c,2) = Tpair.Jaccard(2); % sHRF vs cHRF
        mats(c,3) = Tpair.Jaccard(3); % cHRFderiv vs cHRF
    end

    titles = {'sHRF vs cHRFderiv', 'sHRF vs cHRF', 'cHRFderiv vs cHRF'};

    for j = 1:3
        nexttile;
        imagesc(mats(:,j));
        colormap(parula);
        colorbar;
        yticks(1:numel(contrast_labels));
        yticklabels(contrast_labels);
        xticks([]);
        title(titles{j});
        clim([0 1]);

        for i = 1:numel(contrast_labels)
            text(1, i, sprintf('%.2f', mats(i,j)), ...
                'HorizontalAlignment', 'center', ...
                'Color', 'w', 'FontWeight', 'bold');
        end
    end

    sgtitle('Pairwise Jaccard overlap of minority regions');
end