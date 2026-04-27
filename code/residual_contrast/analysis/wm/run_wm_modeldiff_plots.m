function run_wm_modeldiff_plots(obsFiles, outdir, nPerm, TR)
% run_wm_modeldiff_plots
% Local downstream model-difference plots for WM.
%
% INPUT
%   obsFiles : struct with fields
%              .cHRF
%              .cHRFderiv
%              .sHRF
%   outdir   : folder for saved figures
%   nPerm    : number of sign-flip permutations
%   TR       : repetition time in seconds

if nargin < 2 || isempty(outdir)
    outdir = fullfile(pwd, 'wm_modeldiff_plots');
end
if nargin < 3 || isempty(nPerm)
    nPerm = 5000;
end
if nargin < 4 || isempty(TR)
    TR = 0.72;
end
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

modelNames = {'cHRF', 'cHRFderiv', 'sHRF'};

reqModels = {'cHRF','cHRFderiv','sHRF'};
for i = 1:numel(reqModels)
    assert(isfield(obsFiles, reqModels{i}), 'obsFiles missing field: %s', reqModels{i});
    assert(isfile(obsFiles.(reqModels{i})), 'File not found: %s', obsFiles.(reqModels{i}));
end

D = cell(1,3);
for m = 1:3
    S = load(obsFiles.(reqModels{m}));
    assert(isfield(S, 'TA_subj_oob'),       '%s missing TA_subj_oob',       obsFiles.(reqModels{m}));
    assert(isfield(S, 'PairDiff_subj_oob'), '%s missing PairDiff_subj_oob', obsFiles.(reqModels{m}));
    assert(isfield(S, 'contrast_labels'),   '%s missing contrast_labels',   obsFiles.(reqModels{m}));

    D{m}.TA_subj_oob       = S.TA_subj_oob;        % [S x T x K]
    D{m}.PairDiff_subj_oob = S.PairDiff_subj_oob;  % [S x T x K]
    D{m}.contrast_labels   = cellstr(S.contrast_labels);
end

% check labels match
refLabels = D{1}.contrast_labels(:);
for m = 2:3
    if ~isequal(refLabels, D{m}.contrast_labels(:))
        error('contrast_labels do not match across WM models.');
    end
end

plot_wm_modeldiff_grid(D, modelNames, 'TA_subj_oob', ...
    'Temporal Accuracy Difference', ...
    fullfile(outdir, 'wm_modeldiff_TA.png'), nPerm, TR);

plot_wm_modeldiff_grid(D, modelNames, 'PairDiff_subj_oob', ...
    'Paired Score Difference', ...
    fullfile(outdir, 'wm_modeldiff_PairDiff.png'), nPerm, TR);

end