%% run_WM_region_clustering_pipeline.m
clear; clc;

%% Paths
data_dir = fullfile(pwd, 'data', 'task_residual');
out_dir  = fullfile(data_dir, 'wm_region_clustering_outputs');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

%% Load atlas for parcel names
% Adjust if your atlas is loaded differently
atlas_obj = load_atlas('canlab2018');
region_names = atlas_obj.labels;

%% Load WM residual datasets
% Each file should contain one 4D array:
% [regions x time x subjects x 8]
%
% Expected stimulus order:
% 1=0bk_body, 2=0bk_faces, 3=0bk_places, 4=0bk_tools,
% 5=2bk_body, 6=2bk_faces, 7=2bk_places, 8=2bk_tools

WMsHRF_file      = fullfile(data_dir, 'WMsHRF.mat');
WMcHRFderiv_file = fullfile(data_dir, 'WMcHRFderiv.mat');
WMcHRF_file      = fullfile(data_dir, 'WMcHRF.mat');

WMsHRF      = load_single_var(WMsHRF_file);
WMcHRFderiv = load_single_var(WMcHRFderiv_file);
WMcHRF      = load_single_var(WMcHRF_file);

fprintf('Loaded WMsHRF:      %s\n', mat2str(size(WMsHRF)));
fprintf('Loaded WMcHRFderiv: %s\n', mat2str(size(WMcHRFderiv)));
fprintf('Loaded WMcHRF:      %s\n', mat2str(size(WMcHRF)));

%% Build WM 2bk - 0bk contrasts
[WMdiff_sHRF, contrast_labels]      = build_wm_2b_minus_0b(WMsHRF);
[WMdiff_cHRFderiv, ~]               = build_wm_2b_minus_0b(WMcHRFderiv);
[WMdiff_cHRF, ~]                    = build_wm_2b_minus_0b(WMcHRF);

models = {WMdiff_sHRF, WMdiff_cHRFderiv, WMdiff_cHRF};
model_names = {'sHRF', 'cHRFderiv', 'cHRF'};

%% Run region clustering for all HRF models and WM contrasts
KRange = 2;
NBoot_cluster = 200;     % for region_cluster consensus matrix
NBoot_stab    = 500;     % for formal stability analysis

results_all = struct();
stability_all = struct();

for m = 1:numel(models)
    WMdiff = models{m};
    model_name = model_names{m};

    fprintf('\n=============================\n');
    fprintf('Running model: %s\n', model_name);
    fprintf('=============================\n');

    results_model = cell(numel(contrast_labels), 1);
    stab_model    = cell(numel(contrast_labels), 1);

    for c = 1:numel(contrast_labels)
        contrast_name = contrast_labels{c};
        fprintf('Contrast: %s (2bk - 0bk)\n', contrast_name);

        % [regions x time x subjects] -> [subjects x regions x time]
        X = permute(WMdiff(:, :, :, c), [3 1 2]);

        % Initial region clustering
        res = region_cluster(X, ...
            'ConditionName', sprintf('WM %s %s (2bk-0bk)', model_name, contrast_name), ...
            'KRange', KRange, ...
            'NBoot', NBoot_cluster, ...
            'StandardizeWithinRegion', true, ...
            'RegionNames', region_names, ...
            'MakePlots', true);

        % Formal bootstrap stability analysis
        stab = region_cluster_stability(X, res, ...
            'NBoot', NBoot_stab, ...
            'StandardizeWithinRegion', true, ...
            'MakePlots', true);

        res.model_name = model_name;
        res.contrast_name = contrast_name;

        stab.model_name = model_name;
        stab.contrast_name = contrast_name;

        results_model{c} = res;
        stab_model{c}    = stab;

        % Save per-contrast outputs
        save(fullfile(out_dir, sprintf('WM_%s_%s_region_cluster.mat', model_name, contrast_name)), ...
            'res', 'stab', '-v7.3');

        % Print minority cluster table
        Tminor = get_minority_regions_from_result(res, atlas_obj);
        fprintf('Minority cluster for %s %s:\n', model_name, contrast_name);
        disp(Tminor);
    end

    results_all.(model_name) = results_model;
    stability_all.(model_name) = stab_model;
end

%% Build summary table from initial clustering
Tsummary = summarize_wm_region_clustering_hrfs(results_all);
disp(Tsummary);

writetable(Tsummary, fullfile(out_dir, 'WM_region_clustering_summary.csv'));

%% Build stability summary table
Tstab = summarize_wm_region_stability(stability_all);
disp(Tstab);

writetable(Tstab, fullfile(out_dir, 'WM_region_stability_summary.csv'));

%% Save everything
save(fullfile(out_dir, 'WM_region_clustering_all_results.mat'), ...
    'results_all', 'stability_all', 'Tsummary', 'Tstab', 'contrast_labels', 'model_names', '-v7.3');

fprintf('\nDone. Outputs saved to:\n%s\n', out_dir);

%% ---------- local functions ----------

function X = load_single_var(matfile)
S = load(matfile);
fn = fieldnames(S);
assert(~isempty(fn), 'No variables found in %s', matfile);
X = S.(fn{1});
end

function [WMdiff, contrast_labels] = build_wm_2b_minus_0b(WM)
assert(ndims(WM) == 4, 'WM must be [regions x time x subjects x 8].');
assert(size(WM, 4) == 8, 'WM must have 8 stimulus slots.');

idx_0 = [1 2 3 4];
idx_2 = [5 6 7 8];

WMdiff = WM(:, :, :, idx_2) - WM(:, :, :, idx_0);
contrast_labels = {'body', 'faces', 'places', 'tools'};
end

function Tminor = get_minority_regions_from_result(res, atlas_obj)
counts = accumarray(res.cluster_labels, 1);
[~, minority_cluster] = min(counts);

minor_local_idx = find(res.cluster_labels == minority_cluster);
minor_orig_idx  = res.keep_idx(minor_local_idx);

minor_names = atlas_obj.labels(minor_orig_idx);

Tminor = table( ...
    repmat(minority_cluster, numel(minor_orig_idx), 1), ...
    minor_orig_idx(:), ...
    minor_names(:), ...
    'VariableNames', {'Cluster', 'ParcelIndex', 'ParcelName'});
end

function T = summarize_wm_region_clustering_hrfs(results_all)
model_names = {'sHRF', 'cHRFderiv', 'cHRF'};
contrast_labels = {'body', 'faces', 'places', 'tools'};

rows = {};
for m = 1:numel(model_names)
    model_name = model_names{m};
    res_model = results_all.(model_name);

    for c = 1:numel(contrast_labels)
        res = res_model{c};
        counts = accumarray(res.cluster_labels, 1);
        counts_str = mat2str(counts(:)');

        rows(end+1,:) = { ...
            model_name, ...
            contrast_labels{c}, ...
            res.best_k, ...
            max(res.mean_silhouette), ...
            counts_str ...
            };
    end
end

T = cell2table(rows, 'VariableNames', ...
    {'Model', 'Contrast', 'BestK', 'MeanSilhouette', 'ClusterCounts'});
end

function T = summarize_wm_region_stability(stability_all)
model_names = {'sHRF', 'cHRFderiv', 'cHRF'};
contrast_labels = {'body', 'faces', 'places', 'tools'};

rows = {};

for m = 1:numel(model_names)
    model_name = model_names{m};
    stab_model = stability_all.(model_name);

    for c = 1:numel(contrast_labels)
        stab = stab_model{c};
        ct = stab.cluster_table;

        for r = 1:height(ct)
            rows(end+1,:) = { ...
                model_name, ...
                contrast_labels{c}, ...
                ct.Cluster(r), ...
                ct.Size(r), ...
                ct.WithinConsensus(r), ...
                ct.BetweenConsensus(r), ...
                ct.ConsensusGap(r), ...
                ct.MeanRegionStability(r), ...
                ct.JaccardMean(r), ...
                ct.JaccardSE(r), ...
                stab.overall_ari_mean ...
                };
        end
    end
end

T = cell2table(rows, 'VariableNames', ...
    {'Model','Contrast','Cluster','Size','WithinConsensus','BetweenConsensus', ...
     'ConsensusGap','MeanRegionStability','JaccardMean','JaccardSE','OverallARI'});
end