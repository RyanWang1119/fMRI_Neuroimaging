function TR = detect_wm_tr(repoRoot, cfg)
%DETECT_WM_TR Determine the WM TR from project code or cfg.

if nargin < 2
    cfg = struct();
end

if isfield(cfg, 'TR') && ~isempty(cfg.TR)
    TR = cfg.TR;
    validate_tr(TR);
    return
end

if nargin < 1 || isempty(repoRoot)
    repoRoot = pwd;
end

candidateFiles = { ...
    fullfile(repoRoot, 'code', 'residual_contrast', 'analysis', 'wm', 'plots.m'), ...
    fullfile(repoRoot, 'code', 'residual_contrast', 'analysis', 'wm', 'make_wm_haufe_heatmaps_allcontrasts.m'), ...
    fullfile(repoRoot, 'code', 'residual_contrast', 'mean_accuracy_.m')};

vals = [];
for i = 1:numel(candidateFiles)
    if ~isfile(candidateFiles{i})
        continue
    end
    txt = fileread(candidateFiles{i});
    tok = regexp(txt, 'TR\s*=\s*([0-9]*\.?[0-9]+)\s*;', 'tokens');
    for k = 1:numel(tok)
        vals(end+1) = str2double(tok{k}{1}); %#ok<AGROW>
    end
end

vals = vals(isfinite(vals));
vals = unique(vals);
if isempty(vals)
    error('Could not determine TR from repository code. Supply cfg.TR explicitly.');
end
if numel(vals) > 1
    error('Found conflicting TR values in repository code: %s. Supply cfg.TR explicitly.', mat2str(vals));
end

TR = vals(1);
validate_tr(TR);
end

function validate_tr(TR)
if ~isscalar(TR) || ~isnumeric(TR) || ~isfinite(TR) || TR <= 0
    error('TR must be a finite positive scalar.');
end
end
