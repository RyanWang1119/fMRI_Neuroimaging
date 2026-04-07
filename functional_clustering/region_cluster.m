function results = region_cluster(X, varargin)
% RUN_REGION_CLUSTERING
%
% Region clustering for residual trajectories.
%
% INPUT
%   X : [nSubjects x nRegions x nTime]
%
% OPTIONAL NAME-VALUE PAIRS
%   'ConditionName'           : label for plots/titles
%   'KRange'                  : vector of candidate k values, default 2:6
%   'NBoot'                   : number of bootstrap iterations, default 200
%   'StandardizeWithinRegion' : true/false, z-score each region across time
%   'RegionNames'             : cellstr of region names
%   'MakePlots'               : true/false, default true
%
% OUTPUT
%   results : struct with clustering outputs

% -----------------------------
% Parse inputs
% -----------------------------
p = inputParser;
addParameter(p, 'ConditionName', 'Condition', @(x)ischar(x) || isstring(x));
addParameter(p, 'KRange', 2:6, @isnumeric);
addParameter(p, 'NBoot', 200, @isscalar);
addParameter(p, 'StandardizeWithinRegion', true, @islogical);
addParameter(p, 'RegionNames', {}, @(x)iscell(x) || isstring(x));
addParameter(p, 'MakePlots', true, @islogical);
parse(p, varargin{:});

cond_name = char(p.Results.ConditionName);
k_range   = p.Results.KRange;
nBoot     = p.Results.NBoot;
do_zscore = p.Results.StandardizeWithinRegion;
region_names = p.Results.RegionNames;
make_plots = p.Results.MakePlots;

% -----------------------------
% Check dimensions
% -----------------------------
assert(ndims(X) == 3, 'X must be [subjects x regions x time].');
[nSub, nReg, nTime] = size(X);

if isempty(region_names)
    region_names = arrayfun(@(r) sprintf('Region %d', r), 1:nReg, 'UniformOutput', false);
end

% -----------------------------
% Group-average region x time matrix
% -----------------------------
X_group = squeeze(mean(X, 1, 'omitnan'));   % [nRegions x nTime]

% Optional standardization within each region across time
if do_zscore
    X_use = zscore(X_group, 0, 2);
else
    X_use = X_group;
end

% Handle any all-NaN or zero-variance rows
bad_rows = any(isnan(X_use), 2) | all(abs(X_use) < eps, 2);
if any(bad_rows)
    warning('%d regions removed due to NaN or zero variance.', sum(bad_rows));
    keep_idx = find(~bad_rows);
    X_use = X_use(~bad_rows, :);
    X_group = X_group(~bad_rows, :);
    region_names = region_names(~bad_rows);
else
    keep_idx = 1:nReg;
end

nReg_keep = size(X_use, 1);

% -----------------------------
% Hierarchical clustering
% -----------------------------
D = pdist(X_use, 'correlation');
Z = linkage(D, 'average');

% -----------------------------
% Choose k using silhouette
% -----------------------------
mean_sil = nan(size(k_range));
cluster_labels_all = cell(size(k_range));

for i = 1:numel(k_range)
    k = k_range(i);
    T = cluster(Z, 'maxclust', k);
    cluster_labels_all{i} = T;

    try
        s = silhouette(X_use, T, 'correlation');
        mean_sil(i) = mean(s);
    catch
        % fallback if silhouette fails
        mean_sil(i) = NaN;
    end
end

% pick best k by mean silhouette
[~, best_idx] = max(mean_sil);
best_k = k_range(best_idx);
T_best = cluster_labels_all{best_idx};

% -----------------------------
% Bootstrap co-clustering stability
% -----------------------------
co_clust = zeros(nReg_keep, nReg_keep);
n_valid_boot = 0;

for b = 1:nBoot
    subj_idx = randsample(nSub, nSub, true);
    Xb = squeeze(mean(X(subj_idx, :, :), 1, 'omitnan')); % [nRegions x nTime]

    if do_zscore
        Xb = zscore(Xb, 0, 2);
    end

    Xb = Xb(keep_idx, :);

    bad_b = any(isnan(Xb), 2) | all(abs(Xb) < eps, 2);
    if any(bad_b)
        continue;
    end

    try
        Db = pdist(Xb, 'correlation');
        Zb = linkage(Db, 'average');
        Tb = cluster(Zb, 'maxclust', best_k);

        same_cluster = double(bsxfun(@eq, Tb, Tb'));
        co_clust = co_clust + same_cluster;
        n_valid_boot = n_valid_boot + 1;
    catch
        continue;
    end
end

if n_valid_boot > 0
    co_clust = co_clust / n_valid_boot;
else
    warning('No valid bootstrap samples completed.');
    co_clust = nan(nReg_keep, nReg_keep);
end

% -----------------------------
% Cluster mean trajectories
% -----------------------------
cluster_means = nan(best_k, size(X_group, 2));
cluster_counts = zeros(best_k, 1);

for k = 1:best_k
    idx = (T_best == k);
    cluster_counts(k) = sum(idx);
    cluster_means(k, :) = mean(X_group(idx, :), 1, 'omitnan');
end

% -----------------------------
% Sort regions by cluster for display
% -----------------------------
[~, sort_idx] = sort(T_best);
X_sorted = X_use(sort_idx, :);
T_sorted = T_best(sort_idx);
region_names_sorted = region_names(sort_idx);

% -----------------------------
% Store outputs
% -----------------------------
results = struct();
results.condition_name = cond_name;
results.X_group = X_group;
results.X_use = X_use;
results.keep_idx = keep_idx;
results.bad_rows_removed = bad_rows;
results.distance = 'correlation';
results.linkage = 'average';
results.Z = Z;
results.k_range = k_range;
results.mean_silhouette = mean_sil;
results.best_k = best_k;
results.cluster_labels = T_best;
results.cluster_means = cluster_means;
results.cluster_counts = cluster_counts;
results.co_clust = co_clust;
results.n_valid_boot = n_valid_boot;
results.region_names = region_names;
results.region_names_sorted = region_names_sorted;
results.sort_idx = sort_idx;
results.X_sorted = X_sorted;
results.T_sorted = T_sorted;

% -----------------------------
% Plots
% -----------------------------
if make_plots
    make_region_clustering_plots(results);
end

end


function make_region_clustering_plots(results)

cond_name = results.condition_name;
k_range = results.k_range;
mean_sil = results.mean_silhouette;
Z = results.Z;
best_k = results.best_k;
cluster_means = results.cluster_means;
cluster_counts = results.cluster_counts;
X_sorted = results.X_sorted;
T_sorted = results.T_sorted;
co_clust = results.co_clust;

% 1) Silhouette vs k
figure('Name', [cond_name ' - silhouette'], 'Color', 'w');
plot(k_range, mean_sil, '-o', 'LineWidth', 1.5);
xlabel('Number of clusters (k)');
ylabel('Mean silhouette');
title([cond_name ': choosing k']);
grid on;

% 2) Dendrogram
%figure('Name', [cond_name ' - dendrogram'], 'Color', 'w');
%dendrogram(Z, 0);
%title(sprintf('%s: hierarchical clustering dendrogram', cond_name));
%xlabel('Regions');
%ylabel('Distance');

% 3) Sorted region x time heatmap
%figure('Name', [cond_name ' - sorted trajectories'], 'Color', 'w');
%imagesc(X_sorted);
%xlabel('Time');
%ylabel('Regions (sorted by cluster)');
%title(sprintf('%s: region trajectories sorted by cluster', cond_name));
%colorbar;
%hold on;
%cluster_change = find(diff(T_sorted) ~= 0);
%for i = 1:numel(cluster_change)
%    yline(cluster_change(i) + 0.5, 'k-', 'LineWidth', 1.2);
%end
%hold off;

% 4) Cluster mean trajectories
figure('Name', [cond_name ' - cluster means'], 'Color', 'w');
hold on;
for k = 1:best_k
    plot(cluster_means(k, :), 'LineWidth', 2, ...
        'DisplayName', sprintf('Cluster %d (n=%d)', k, cluster_counts(k)));
end
xlabel('Time');
ylabel('Mean residual');
title(sprintf('%s: cluster-average residual trajectories', cond_name));
legend('Location', 'best');
grid on;
hold off;

% 5) Bootstrap co-clustering matrix
%figure('Name', [cond_name ' - coclustering'], 'Color', 'w');
%imagesc(co_clust);
%axis square;
%xlabel('Regions');
%ylabel('Regions');
%title(sprintf('%s: bootstrap co-clustering stability', cond_name));
%colorbar;

end

