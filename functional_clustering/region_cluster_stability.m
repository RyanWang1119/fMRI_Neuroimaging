function stab = region_cluster_stability(X, res, varargin)
% REGION_CLUSTER_STABILITY
%
% Formal bootstrap stability analysis for region clustering.
%
% INPUT
%   X   : [subjects x regions x time]
%   res : output struct from region_cluster(X, ...)
%
% OPTIONAL NAME-VALUE PAIRS
%   'NBoot'                   : number of bootstrap iterations, default 500
%   'StandardizeWithinRegion' : true/false, default true
%   'MakePlots'               : true/false, default true
%
% OUTPUT
%   stab : struct containing:
%       .overall_ari_mean
%       .overall_ari_se
%       .cluster_table
%       .consensus_sorted
%       .region_stability
%       .boot_ari
%       .boot_jaccard
%
% Notes:
% - Uses the reference clustering already stored in res.cluster_labels
% - Re-clusters bootstrap subject averages with the same best_k
% - ARI is label-invariant
% - Cluster-wise Jaccard is computed after optimal label matching

p = inputParser;
addParameter(p, 'NBoot', 500, @isscalar);
addParameter(p, 'StandardizeWithinRegion', true, @islogical);
addParameter(p, 'MakePlots', true, @islogical);
parse(p, varargin{:});

nBoot = p.Results.NBoot;
do_zscore = p.Results.StandardizeWithinRegion;
make_plots = p.Results.MakePlots;

assert(ndims(X) == 3, 'X must be [subjects x regions x time].');

[nSub, ~, ~] = size(X);

% Reference solution
T_ref = res.cluster_labels(:);
keep_idx = res.keep_idx(:);
k = res.best_k;
cond_name = res.condition_name;

nReg_keep = numel(T_ref);

% Bootstrap outputs
co_clust = zeros(nReg_keep, nReg_keep);
boot_ari = nan(nBoot, 1);
boot_jaccard = nan(nBoot, k);

n_valid = 0;

for b = 1:nBoot
    subj_idx = randsample(nSub, nSub, true);

    % group-average region x time
    Xb = squeeze(mean(X(subj_idx, :, :), 1, 'omitnan'));   % [regions x time]

    if do_zscore
        Xb = zscore(Xb, 0, 2);
    end

    Xb = Xb(keep_idx, :);

    bad_rows = any(isnan(Xb), 2) | all(abs(Xb) < eps, 2);
    if any(bad_rows)
        continue;
    end

    try
        Db = pdist(Xb, 'correlation');
        Zb = linkage(Db, 'average');
        Tb = cluster(Zb, 'maxclust', k);

        n_valid = n_valid + 1;

        % Consensus matrix
        co_clust = co_clust + double(bsxfun(@eq, Tb, Tb'));

        % Overall stability
        boot_ari(n_valid) = adjusted_rand_index(T_ref, Tb);

        % Cluster-wise matched Jaccard
        boot_jaccard(n_valid, :) = matched_cluster_jaccard(T_ref, Tb, k);

    catch
        continue;
    end
end

% Trim to valid runs
boot_ari = boot_ari(1:n_valid);
boot_jaccard = boot_jaccard(1:n_valid, :);

if n_valid == 0
    error('No valid bootstrap runs completed.');
end

co_clust = co_clust ./ n_valid;

% Sort consensus by reference cluster
[~, sort_idx] = sort(T_ref);
consensus_sorted = co_clust(sort_idx, sort_idx);
T_sorted = T_ref(sort_idx);

% Region-level stability:
% mean consensus with members of its own reference cluster
region_stability = nan(nReg_keep, 1);
for i = 1:nReg_keep
    same = find(T_ref == T_ref(i));
    same(same == i) = [];
    if ~isempty(same)
        region_stability(i) = mean(co_clust(i, same));
    end
end

% Cluster-level summary table
cluster_id = (1:k)';
cluster_size = accumarray(T_ref, 1, [k 1]);
within_consensus = nan(k,1);
between_consensus = nan(k,1);
stability_gap = nan(k,1);
mean_region_stability = nan(k,1);
jaccard_mean = nan(k,1);
jaccard_se = nan(k,1);

for c = 1:k
    idx = find(T_ref == c);
    not_idx = find(T_ref ~= c);

    if numel(idx) >= 2
        A = co_clust(idx, idx);
        mask = triu(true(numel(idx)), 1);
        within_consensus(c) = mean(A(mask), 'omitnan');
    else
        within_consensus(c) = NaN;
    end

    if ~isempty(not_idx)
        between_consensus(c) = mean(co_clust(idx, not_idx), 'all', 'omitnan');
    end

    stability_gap(c) = within_consensus(c) - between_consensus(c);
    mean_region_stability(c) = mean(region_stability(idx), 'omitnan');

    jaccard_mean(c) = mean(boot_jaccard(:, c), 'omitnan');
    jaccard_se(c)   = std(boot_jaccard(:, c), 0, 1, 'omitnan') ./ sqrt(sum(~isnan(boot_jaccard(:, c))));
end

cluster_table = table( ...
    cluster_id, ...
    cluster_size, ...
    within_consensus, ...
    between_consensus, ...
    stability_gap, ...
    mean_region_stability, ...
    jaccard_mean, ...
    jaccard_se, ...
    'VariableNames', { ...
    'Cluster', ...
    'Size', ...
    'WithinConsensus', ...
    'BetweenConsensus', ...
    'ConsensusGap', ...
    'MeanRegionStability', ...
    'JaccardMean', ...
    'JaccardSE'});

stab = struct();
stab.condition_name = cond_name;
stab.n_valid_boot = n_valid;
stab.boot_ari = boot_ari;
stab.overall_ari_mean = mean(boot_ari, 'omitnan');
stab.overall_ari_se = std(boot_ari, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(boot_ari)));
stab.boot_jaccard = boot_jaccard;
stab.cluster_table = cluster_table;
stab.consensus = co_clust;
stab.consensus_sorted = consensus_sorted;
stab.sort_idx = sort_idx;
stab.T_sorted = T_sorted;
stab.region_stability = region_stability;

if make_plots
    make_region_cluster_stability_plots(stab);
end

end


% ============================================================
% Plotting
% ============================================================
function make_region_cluster_stability_plots(stab)

cond_name = stab.condition_name;
cluster_table = stab.cluster_table;
consensus_sorted = stab.consensus_sorted;
T_sorted = stab.T_sorted;
boot_ari = stab.boot_ari;

% 1) Sorted consensus matrix
figure('Name', [cond_name ' - stability consensus'], 'Color', 'w');
imagesc(consensus_sorted);
axis square;
xlabel('Regions (sorted by reference cluster)');
ylabel('Regions (sorted by reference cluster)');
title(sprintf('%s: bootstrap consensus matrix', cond_name));
colorbar;
hold on;
cluster_change = find(diff(T_sorted) ~= 0);
for i = 1:numel(cluster_change)
    xline(cluster_change(i) + 0.5, 'k-', 'LineWidth', 1.2);
    yline(cluster_change(i) + 0.5, 'k-', 'LineWidth', 1.2);
end
hold off;

% 2) Cluster-wise Jaccard
figure('Name', [cond_name ' - cluster jaccard'], 'Color', 'w');
bar(cluster_table.Cluster, cluster_table.JaccardMean);
hold on;
errorbar(cluster_table.Cluster, cluster_table.JaccardMean, cluster_table.JaccardSE, ...
    'k.', 'LineWidth', 1.2);
xlabel('Cluster');
ylabel('Bootstrap Jaccard');
title(sprintf('%s: cluster-wise bootstrap Jaccard', cond_name));
grid on;
hold off;

% 3) ARI histogram
figure('Name', [cond_name ' - ARI'], 'Color', 'w');
histogram(boot_ari);
xlabel('Adjusted Rand Index');
ylabel('Bootstrap count');
title(sprintf('%s: overall clustering stability (ARI)', cond_name));
grid on;

end


% ============================================================
% Adjusted Rand Index
% ============================================================
function ari = adjusted_rand_index(labels1, labels2)

labels1 = labels1(:);
labels2 = labels2(:);

assert(numel(labels1) == numel(labels2), 'Label vectors must have same length.');

u1 = unique(labels1);
u2 = unique(labels2);

n1 = numel(u1);
n2 = numel(u2);

cont = zeros(n1, n2);

for i = 1:n1
    for j = 1:n2
        cont(i,j) = sum(labels1 == u1(i) & labels2 == u2(j));
    end
end

nij2 = sum(cont(:) .* (cont(:) - 1) / 2);
ai2  = sum(sum(cont, 2) .* (sum(cont, 2) - 1) / 2);
bj2  = sum(sum(cont, 1) .* (sum(cont, 1) - 1) / 2);
n = sum(cont(:));
n2tot = n * (n - 1) / 2;

expected = (ai2 * bj2) / n2tot;
max_index = 0.5 * (ai2 + bj2);

denom = max_index - expected;
if abs(denom) < eps
    ari = NaN;
else
    ari = (nij2 - expected) / denom;
end

end


% ============================================================
% Cluster-wise Jaccard with optimal label matching
% ============================================================
function jbest = matched_cluster_jaccard(T_ref, T_boot, k)

J = zeros(k, k);

for c = 1:k
    A = (T_ref == c);
    for d = 1:k
        B = (T_boot == d);
        inter = sum(A & B);
        union = sum(A | B);
        if union == 0
            J(c,d) = NaN;
        else
            J(c,d) = inter / union;
        end
    end
end

% Brute-force optimal assignment; fine because k is small (2-6)
perm_list = perms(1:k);
scores = nan(size(perm_list,1), 1);

for p = 1:size(perm_list,1)
    idx = sub2ind([k k], 1:k, perm_list(p,:));
    scores(p) = sum(J(idx), 'omitnan');
end

[~, best_p] = max(scores);
best_perm = perm_list(best_p, :);

jbest = nan(1, k);
for c = 1:k
    jbest(c) = J(c, best_perm(c));
end

end