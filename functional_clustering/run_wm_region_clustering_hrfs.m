function results_all = run_wm_region_clustering_hrfs(WMsHRF, WMcHRFderiv, WMcHRF, varargin)
% RUN_WM_REGION_CLUSTERING_HRFS
%
% INPUTS
%   WMsHRF      : [regions x time x subjects x 8]
%   WMcHRFderiv : [regions x time x subjects x 8]
%   WMcHRF      : [regions x time x subjects x 8]
%
%   'RegionNames' : cell array of region names
%   'KRange'      : candidate k values, default 2:6
%   'NBoot'       : bootstrap iterations, default 200
%   'MakePlots'   : true/false, default true
%
% OUTPUT
%   results_all.(model_name){contrast_idx}

p = inputParser;
addParameter(p, 'RegionNames', {}, @(x) iscell(x) || isstring(x));
addParameter(p, 'KRange', 2:6, @isnumeric);
addParameter(p, 'NBoot', 200, @isscalar);
addParameter(p, 'MakePlots', true, @islogical);
parse(p, varargin{:});

region_names = p.Results.RegionNames;
k_range = p.Results.KRange;
nBoot = p.Results.NBoot;
make_plots = p.Results.MakePlots;

models = {WMsHRF, WMcHRFderiv, WMcHRF};
model_names = {'sHRF', 'cHRFderiv', 'cHRF'};

results_all = struct();

for m = 1:numel(models)
    WM = models{m};

    [WMdiff, contrast_labels] = build_wm_contrast(WM);

    results_model = cell(numel(contrast_labels), 1);

    fprintf('\n=== Running %s ===\n', model_names{m});

    for c = 1:numel(contrast_labels)
        fprintf('Contrast: %s (2bk - 0bk)\n', contrast_labels{c});

        X = permute(WMdiff(:, :, :, c), [3 1 2]);

        res = region_cluster(X, ...
            'ConditionName', sprintf('WM %s: %s (2bk-0bk)', model_names{m}, contrast_labels{c}), ...
            'KRange', k_range, ...
            'NBoot', nBoot, ...
            'StandardizeWithinRegion', true, ...
            'RegionNames', region_names, ...
            'MakePlots', make_plots);

        res.model_name = model_names{m};
        res.contrast_name = contrast_labels{c};
        res.original_size = size(WMdiff(:, :, :, c));
        res.reordered_size = size(X);

        results_model{c} = res;
    end

    results_all.(model_names{m}) = results_model;
end

results_all.contrast_labels = contrast_labels;
results_all.model_names = model_names;

end