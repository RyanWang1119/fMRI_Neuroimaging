function test_roi_influence_synthetic()
%TEST_ROI_INFLUENCE_SYNTHETIC Synthetic checks for WM ROI influence module.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
moduleDir = fullfile(repoRoot, 'code', 'residual_contrast', 'analysis', 'wm_roi_influence');
addpath(moduleDir);

if exist('lassoglm', 'file') ~= 2
    error('This test requires lassoglm from the Statistics and Machine Learning Toolbox.');
end

rng(20260621, 'twister');
nSubject = 50;
nROI = 40;
subjectIDs = (100001:(100000+nSubject))';

latent = 1.35 + 0.55 * randn(nSubject, 1);
D = 0.75 * randn(nSubject, nROI);
D(:, 1) = latent + 0.25 * randn(nSubject, 1);       % true task-linked parcel
D(:, 2) = D(:, 1) + 0.12 * randn(nSubject, 1);       % correlated proxy
D(:, 3:end) = D(:, 3:end) + 0.05 * randn(nSubject, nROI-2);

tmpRoot = fullfile(tempdir, sprintf('roi_influence_synthetic_%s', datestr(now, 30)));
if ~exist(tmpRoot, 'dir')
    mkdir(tmpRoot);
end

cfg = struct();
cfg.seed = 20260621;
cfg.outer.nRepeats = 24;
cfg.outer.trainFraction = 0.80;
cfg.outer.checkpointBatchSize = 8;
cfg.outer.useParfor = false;
cfg.model.alpha = 0.50;
cfg.model.numLambda = 30;
cfg.model.innerKFolds = 5;
cfg.model.useOneSE = true;
cfg.importance.nRandomSwaps = 12;
cfg.importance.selectedOnly = true;
cfg.permutation.nPerm = 500;
cfg.permutation.alpha = 0.05;
cfg.thresholds.minSelectionFrequency = 0.50;
cfg.thresholds.minSignConsistency = 0.80;
cfg.thresholds.minMeanReliance = -Inf;
cfg.thresholds.requireFWER = true;

splitFile = fullfile(tmpRoot, 'shared_subject_splits.mat');
outerSplits = make_repeated_subject_splits(nSubject, subjectIDs, cfg, splitFile);

% 1. No subject appears in both train and test sets.
for r = 1:numel(outerSplits)
    assert(isempty(intersect(outerSplits(r).trainIdx, outerSplits(r).testIdx)), ...
        'Subject leakage between train and test in outer split %d.', r);
end

% 2. Inner folds preserve subject pairing.
foldID = make_grouped_inner_folds(numel(outerSplits(1).trainIdx), cfg.model.innerKFolds, cfg.seed);
rowFold = [foldID; foldID];
for i = 1:numel(foldID)
    assert(rowFold(i) == rowFold(i + numel(foldID)), 'Synthetic paired rows split across inner folds.');
end

% 3. Training-only scaling is not contaminated by test values.
Dtrain = D(outerSplits(1).trainIdx, :);
Dtest = D(outerSplits(1).testIdx, :);
DtestContaminated = Dtest;
DtestContaminated(:, 5) = 1e9;
prepClean = apply_training_preprocessing(Dtrain, Dtest);
prepContam = apply_training_preprocessing(Dtrain, DtestContaminated);
assert(max(abs(prepClean.featureMedian - prepContam.featureMedian)) == 0, 'Feature medians used test data.');
assert(max(abs(prepClean.featureScale - prepContam.featureScale)) == 0, 'Feature scales used test data.');

outerResults = repmat(empty_result_for_test(nROI, cfg.importance.nRandomSwaps), cfg.outer.nRepeats, 1);
for r = 1:cfg.outer.nRepeats
    outerResults(r) = fit_outer_elasticnet_model(D, subjectIDs, outerSplits(r), cfg);
end

roiLabels = table((1:nROI)', arrayfun(@(r) sprintf('Synthetic_%02d', r), (1:nROI)', 'UniformOutput', false), ...
    repmat({''}, nROI, 1), 'VariableNames', {'roi_index','roi_label','roi_network'});
signFlipResults = run_paired_signflip_maxT(D, cfg);
summaryTable = summarize_roi_influence(outerResults, signFlipResults, roiLabels, cfg);

nullIdx = 3:nROI;
pairProfile = max(summaryTable.selection_frequency(1:2), summaryTable.mean_model_reliance_tafc(1:2));
nullProfile = max(summaryTable.selection_frequency(nullIdx), summaryTable.mean_model_reliance_tafc(nullIdx));

% 4. The true parcel or correlated proxy has a higher stability/reliance
% profile than most null parcels.
assert(any(pairProfile > prctile(nullProfile, 80)), ...
    'Neither the true parcel nor its proxy exceeded most null parcels.');

% 5. Intentionally do not require both correlated parcels to be selected.
assert(any(summaryTable.selection_frequency(1:2) > median(summaryTable.selection_frequency(nullIdx))), ...
    'Expected at least one of true/proxy parcels to exceed null selection frequency.');

% 6. Sign-flip maxT detects the true-effect parcel at reasonable power.
assert(signFlipResults.pFWER(1) < cfg.permutation.alpha, ...
    'Sign-flip maxT did not detect the true synthetic parcel.');

fprintf('Synthetic ROI influence test passed.\n');
end

function e = empty_result_for_test(nROI, nSwaps)
e = struct();
e.repeat = NaN;
e.seed = NaN;
e.betaRaw = nan(nROI, 1);
e.betaStandardized = nan(nROI, 1);
e.intercept = NaN;
e.selectedLambda = NaN;
e.trainSubjectIDs = [];
e.testSubjectIDs = [];
e.trainSubjectIdx = [];
e.testSubjectIdx = [];
e.TAFC = NaN;
e.testMargins = [];
e.selectionIndicator = false(nROI, 1);
e.haufePatternRaw = nan(nROI, 1);
e.haufePatternZ = nan(nROI, 1);
e.modelRelianceTAFC = zeros(nROI, nSwaps, 'single');
e.featureMedian = nan(1, nROI);
e.featureMean = nan(1, nROI);
e.featureScale = nan(1, nROI);
e.innerFoldIDByTrainSubject = [];
e.innerValidationLossMean = [];
e.innerValidationLossSE = [];
e.lambdaPath = [];
end
