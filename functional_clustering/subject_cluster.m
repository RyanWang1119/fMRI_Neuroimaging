function results = subject_cluster(X, varargin)
% SUBJECT_CLUSTER
%
% Cluster subjects based on region x time residual signatures.
%
% INPUT
%   X : [nSubjects x nRegions x nTime]
%
% OPTIONAL NAME-VALUE PAIRS
%   'ConditionName'              : label for plots/titles
%   'KRange'                     : vector of candidate k values, default 2:6
%   'NBoot'                      : number of bootstrap iterations, default 200
%   'FeatureMode'                : 'flatten' or 'mean_time', default 'flatten'
%   'StandardizeWithinSubject'   : true/false, default true
%   'DoPCA'                      : true/false, default true
%   'NumPCs'                     : fixed number of PCs to keep, default []
%   'VarianceToKeep'             : cumulative % variance if NumPCs empty, default 90
%   'MakePlots'                  : true/false, default true
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
addParameter(p, 'FeatureMode', 'flatten', @(x)ischar(x) || isstring(x));
addParameter(p, 'StandardizeWithinSubject', true, @islogical);
addParameter(p, 'DoPCA', true, @islogical);
addParameter(p, 'NumPCs', [], @(x)isempty(x) || isscalar(x));
addParameter(p, 'VarianceToKeep', 90, @isscalar);
addParameter(p, 'MakePlots', true, @islogical);
parse(p, varargin{:});

cond_name   = char(p.Results.ConditionName);
k_range     = p.Results.KRange;
nBoot       = p.Results.NBoot;
feature_mode = lower(char(p.Results.FeatureMode));
do_row_z    = p.Results.StandardizeWithinSubject;
do_pca      = p.Results.DoPCA;
numPCs      = p.Results.NumPCs;
var_keep    = p.Results.VarianceToKeep;
make_plots  = p.Results.MakePlots;

% -----------------------------
% Check input
% -----------------------------
assert(ndims(X) == 3, 'X must be [subjects x regions x time].');
[nSub, nReg, nTime] = size(X);

% -----------------------------
% Build subject feature matrix
% -----------------------------
[X_feat, X_meantraj] = build_subject_features(X, feature_mode);

% Row-standardize each subject if desired
if do_row_z
    X_feat = safe_row_zscore(X_feat);
end

% Remove bad columns (constant or NaN across subjects)
bad_cols = any(isnan(X_feat), 1) | (std(X_feat, 0, 1) < eps);
X_feat = X_feat(:, ~bad_cols);

% Remove bad subjects if needed
bad_rows = any(isnan(X_feat), 2) | all(abs(X_feat) < eps, 2);
if any(bad_rows)
    warning('%d subjects removed due to NaN or zero variance.', sum(bad_rows));
end
keep_idx = find(~bad_rows);

X_feat = X_feat(keep_idx, :);
X_meantraj = X_meantraj(keep_idx, :);
X_base = X(keep_idx, :, :);

nSub_keep = size(X_feat, 1);

% -----------------------------
% PCA / clustering representation
% -----------------------------
if do_pca
    [coeff, score, latent, ~, explained, mu] = pca(X_feat, 'Economy', true);

    if isempty(numPCs)
        nPC = find(cumsum(explained) >= var_keep, 1, 'first');
        if isempty(nPC)
            nPC = min(size(score,2), 10);
        end
        nPC = max(2, min(nPC, size(score,2)));
    else
        nPC = min(numPCs, size(score,2));
    end

    X_cluster = score(:, 1:nPC);
    distance_name = 'euclidean';
    linkage_name  = 'ward';
else
    coeff = [];
    score = [];
    latent = [];
    explained = [];
    mu = [];
    nPC = [];
    X_cluster = X_feat;
    distance_name = 'correlation';
    linkage_name  = 'average';
end

% -----------------------------
% Hierarchical clustering
% -----------------------------
D = pdist(X_cluster, distance_name);
Z = linkage(D, linkage_name);

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
        s = silhouette(X_cluster, T, distance_name);
        mean_sil(i) = mean(s);
    catch
        mean_sil(i) = NaN;
    end
end

[~, best_idx] = max(mean_sil);
best_k = k_range(best_idx);
T_best = cluster_labels_all{best_idx};

% -----------------------------
% Bootstrap co-clustering stability
% Resample regions with replacement
% -----------------------------
co_clust = zeros(nSub_keep, nSub_keep);
n_valid_boot = 0;

for b = 1:nBoot
    reg_idx = randsample(nReg, nReg, true);
    Xb = X_base(:, reg_idx, :);

    [Xb_feat, ~] = build_subject_features(Xb, feature_mode);

    if do_row_z
        Xb_feat = safe_row_zscore(Xb_feat);
    end

    bad_cols_b = any(isnan(Xb_feat), 1) | (std(Xb_feat, 0, 1) < eps);
    Xb_feat = Xb_feat(:, ~bad_cols_b);

    bad_rows_b = any(isnan(Xb_feat), 2) | all(abs(Xb_feat) < eps, 2);
    if any(bad_rows_b)
        continue;
    end

    try
        if do_pca
            [~, score_b, ~, ~, explained_b] = pca(Xb_feat, 'Economy', true);

            if isempty(numPCs)
                nPC_b = find(cumsum(explained_b) >= var_keep, 1, 'first');
                if isempty(nPC_b)
                    nPC_b = min(size(score_b,2), 10);
                end
                nPC_b = max(2, min(nPC_b, size(score_b,2)));
            else
                nPC_b = min(numPCs, size(score_b,2));
            end

            Xb_cluster = score_b(:, 1:nPC_b);
            Db = pdist(Xb_cluster, 'euclidean');
            Zb = linkage(Db, 'ward');
        else
            Db = pdist(Xb_feat, 'correlation');
            Zb = linkage(Db, 'average');
        end

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
    co_clust = nan(nSub_keep, nSub_keep);
end

% -----------------------------
% Cluster mean trajectories (subject-average over regions)
% -----------------------------
cluster_meantraj = nan(best_k, size(X_meantraj, 2));
cluster_counts = zeros(best_k, 1);

for k = 1:best_k
    idx = (T_best == k);
    cluster_counts(k) = sum(idx);
    cluster_meantraj(k, :) = mean(X_meantraj(idx, :), 1, 'omitnan');
end

% -----------------------------
% Sort subjects by cluster for display
% -----------------------------
[~, sort_idx] = sort(T_best);
T_sorted = T_best(sort_idx);
X_meantraj_sorted = X_meantraj(sort_idx, :);

if do_pca && size(X_cluster,2) >= 2
    X_scatter = X_cluster(:,1:2);
else
    X_scatter = [];
end

% -----------------------------
% Store outputs
% -----------------------------
results = struct();
results.condition_name = cond_name;
results.feature_mode = feature_mode;
results.n_subjects = nSub;
results.n_regions = nReg;
results.n_time = nTime;
results.keep_idx = keep_idx;
results.bad_rows_removed = bad_rows;
results.distance = distance_name;
results.linkage = linkage_name;
results.Z = Z;
results.k_range = k_range;
results.mean_silhouette = mean_sil;
results.best_k = best_k;
results.cluster_labels = T_best;
results.cluster_counts = cluster_counts;
results.cluster_meantraj = cluster_meantraj;
results.co_clust = co_clust;
results.n_valid_boot = n_valid_boot;
results.X_feat = X_feat;
results.X_cluster = X_cluster;
results.X_meantraj = X_meantraj;
results.X_meantraj_sorted = X_meantraj_sorted;
results.T_sorted = T_sorted;
results.sort_idx = sort_idx;

if do_pca
    results.pca_coeff = coeff;
    results.pca_score = score;
    results.pca_latent = latent;
    results.pca_explained = explained;
    results.pca_mu = mu;
    results.nPC = nPC;
else
    results.pca_coeff = [];
    results.pca_score = [];
    results.pca_latent = [];
    results.pca_explained = [];
    results.pca_mu = [];
    results.nPC = [];
end

% -----------------------------
% Plots
% -----------------------------
if make_plots
    make_subject_clustering_plots(results, X_scatter);
end

end


% ============================================================
% Helper: build subject features
% ============================================================
function [X_feat, X_meantraj] = build_subject_features(X, feature_mode)
% X is [subjects x regions x time]
[nSub, nReg, nTime] = size(X);

switch feature_mode
    case 'flatten'
        % subject x (region*time)
        X_feat = reshape(X, nSub, nReg * nTime);

    case 'mean_time'
        % average over regions -> subject x time
        X_feat = squeeze(mean(X, 2, 'omitnan'));
        if isvector(X_feat)
            X_feat = X_feat(:);
        end

    otherwise
        error('Unknown FeatureMode: %s. Use ''flatten'' or ''mean_time''.', feature_mode);
end

% always keep this for interpretation
X_meantraj = squeeze(mean(X, 2, 'omitnan'));
if isvector(X_meantraj)
    X_meantraj = X_meantraj(:);
end
end


% ============================================================
% Helper: safe row z-score
% ============================================================
function Xz = safe_row_zscore(X)
mu = mean(X, 2, 'omitnan');
sd = std(X, 0, 2, 'omitnan');
sd(sd < eps | isnan(sd)) = 1;

Xz = (X - mu) ./ sd;
Xz(~isfinite(Xz)) = 0;
end


% ============================================================
% Helper: plots
% ============================================================
function make_subject_clustering_plots(results, X_scatter)

cond_name = results.condition_name;
k_range = results.k_range;
mean_sil = results.mean_silhouette;
Z = results.Z;
best_k = results.best_k;
cluster_meantraj = results.cluster_meantraj;
cluster_counts = results.cluster_counts;
X_meantraj_sorted = results.X_meantraj_sorted;
T_sorted = results.T_sorted;
co_clust = results.co_clust;
X_cluster = results.X_cluster;
cluster_labels = results.cluster_labels;

% 1) Silhouette vs k
figure('Name', [cond_name ' - subject silhouette'], 'Color', 'w');
plot(k_range, mean_sil, '-o', 'LineWidth', 1.5);
xlabel('Number of clusters (k)');
ylabel('Mean silhouette');
title([cond_name ': subject clustering - choosing k']);
grid on;

% 2) Dendrogram
figure('Name', [cond_name ' - subject dendrogram'], 'Color', 'w');
dendrogram(Z, 0);
title(sprintf('%s: subject clustering dendrogram', cond_name));
xlabel('Subjects');
ylabel('Distance');

% 3) PCA scatter if available
if ~isempty(X_scatter)
    figure('Name', [cond_name ' - subject PCA scatter'], 'Color', 'w');
    gscatter(X_scatter(:,1), X_scatter(:,2), cluster_labels);
    xlabel('PC 1');
    ylabel('PC 2');
    title(sprintf('%s: subject clusters in PCA space', cond_name));
    grid on;
end

% 4) Sorted subject-average trajectories
figure('Name', [cond_name ' - sorted subject mean trajectories'], 'Color', 'w');
imagesc(X_meantraj_sorted);
xlabel('Time');
ylabel('Subjects (sorted by cluster)');
title(sprintf('%s: subject mean trajectories sorted by cluster', cond_name));
colorbar;
hold on;
cluster_change = find(diff(T_sorted) ~= 0);
for i = 1:numel(cluster_change)
    yline(cluster_change(i) + 0.5, 'k-', 'LineWidth', 1.2);
end
hold off;

% 5) Cluster-average subject mean trajectories
figure('Name', [cond_name ' - subject cluster means'], 'Color', 'w');
hold on;
for k = 1:best_k
    plot(cluster_meantraj(k, :), 'LineWidth', 2, ...
        'DisplayName', sprintf('Cluster %d (n=%d)', k, cluster_counts(k)));
end
xlabel('Time');
ylabel('Mean residual across regions');
title(sprintf('%s: subject cluster-average trajectories', cond_name));
legend('Location', 'best');
grid on;
hold off;

% 6) Bootstrap co-clustering matrix
figure('Name', [cond_name ' - subject coclustering'], 'Color', 'w');
imagesc(co_clust);
axis square;
xlabel('Subjects');
ylabel('Subjects');
title(sprintf('%s: subject bootstrap co-clustering stability', cond_name));
colorbar;

% 7) Optional: heatmap of clustering features if low-dimensional
if size(X_cluster, 2) <= 30
    [~, sort_idx] = sort(cluster_labels);
    figure('Name', [cond_name ' - subject feature heatmap'], 'Color', 'w');
    imagesc(X_cluster(sort_idx, :));
    xlabel('Clustering features');
    ylabel('Subjects (sorted by cluster)');
    title(sprintf('%s: clustering representation', cond_name));
    colorbar;
end

end