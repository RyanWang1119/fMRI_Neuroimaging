function run_language_modeldiff_plots(baseDir, outdir, nPerm)

if nargin < 1 || isempty(baseDir)
    baseDir = 'C:\Users\Ryan_\Documents\my_Git\fMRI\data\observed_residual_results';
end
if nargin < 2 || isempty(outdir)
    outdir = fullfile(baseDir, 'language_simple_analysis');
end
if nargin < 3 || isempty(nPerm)
    nPerm = 5000;
end
if ~exist(outdir, 'dir'), mkdir(outdir); end

files = { ...
    fullfile(baseDir,'observed_language_simple_chrf.mat'), ...
    fullfile(baseDir,'observed_language_simple_chrfderiv.mat'), ...
    fullfile(baseDir,'observed_language_simple_shrf.mat')};

modelNames = {'cHRF','cHRFderiv','sHRF'};

D = cell(1,3);
for m = 1:3
    S = load(files{m});
    D{m}.TA_subj_oob       = S.TA_subj_oob;
    D{m}.PairDiff_subj_oob = S.PairDiff_subj_oob;
    D{m}.contrast_labels   = cellstr(S.contrast_labels);
end

plot_language_modeldiff_3contrasts(D, modelNames, 'TA_subj_oob', ...
    'Temporal Accuracy Difference', ...
    'Model Differences in Temporal Accuracy: Phase-Matched Story vs. Math Contrasts', ...
    fullfile(outdir, 'language_modeldiff_TA_3contrasts.png'), nPerm);

plot_language_modeldiff_3contrasts(D, modelNames, 'PairDiff_subj_oob', ...
    'Paired Score Difference', ...
    'Model Differences in Paired Score Difference: Phase-Matched Story vs. Math Contrasts', ...
    fullfile(outdir, 'language_modeldiff_PairDiff_3contrasts.png'), nPerm);

end
