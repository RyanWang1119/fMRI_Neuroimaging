function result = fit_outer_elasticnet_model(D, subjectIDs, outerSplit, cfg)
%FIT_OUTER_ELASTICNET_MODEL Fit one outer split and evaluate held-out TAFC.

if exist('lassoglm', 'file') ~= 2
    error('lassoglm is required for elastic-net logistic regression.');
end

D = double(D);
subjectIDs = subjectIDs(:);
[nSubject, nROI] = size(D);

trainIdx = outerSplit.trainIdx(:);
testIdx = outerSplit.testIdx(:);
if ~isempty(intersect(trainIdx, testIdx))
    error('Outer repeat %d has overlapping train/test subjects.', outerSplit.repeat);
end
if any(trainIdx < 1 | trainIdx > nSubject) || any(testIdx < 1 | testIdx > nSubject)
    error('Outer repeat %d has subject indices outside range.', outerSplit.repeat);
end

Dtrain = D(trainIdx, :);
Dtest = D(testIdx, :);

foldSeed = outerSplit.seed + 17;
foldID = make_grouped_inner_folds(numel(trainIdx), cfg.model.innerKFolds, foldSeed);
tune = tune_elasticnet_logistic_groupcv(Dtrain, cfg, foldID);

prep = apply_training_preprocessing(Dtrain, Dtest);
ypairTrain = [ones(numel(trainIdx), 1); zeros(numel(trainIdx), 1)];

[B, fitInfo] = lassoglm(prep.XpairTrainZ, ypairTrain, 'binomial', ...
    'Alpha', cfg.model.alpha, ...
    'Lambda', tune.selectedLambda, ...
    'Standardize', false);

betaStandardized = B(:, 1);
intercept = fitInfo.Intercept(1);

% betaRaw is in original paired-row feature units. Because D margins compare
% +D/2 against -D/2, held-out margins in original D units are D * betaRaw.
betaRaw = betaStandardized ./ prep.featureScale(:);

testMargins = prep.DapplyZ * betaStandardized;
TAFC = mean(testMargins > 0, 'omitnan');

selectionIndicator = betaStandardized ~= 0;
[haufePatternZ, haufePatternRaw] = compute_haufe_pattern( ...
    prep.XpairTrainZ, betaStandardized, prep.featureScale);

importanceSeed = outerSplit.seed + 313;
[modelRelianceTAFC, ~] = compute_heldout_model_reliance( ...
    prep.DapplyZ, betaStandardized, selectionIndicator, cfg, importanceSeed);

result = struct();
result.repeat = outerSplit.repeat;
result.seed = outerSplit.seed;
result.betaRaw = betaRaw(:);
result.betaStandardized = betaStandardized(:);
result.intercept = intercept;
result.selectedLambda = tune.selectedLambda;
result.trainSubjectIDs = subjectIDs(trainIdx);
result.testSubjectIDs = subjectIDs(testIdx);
result.trainSubjectIdx = trainIdx;
result.testSubjectIdx = testIdx;
result.TAFC = TAFC;
result.testMargins = testMargins(:);
result.selectionIndicator = selectionIndicator(:);
result.haufePatternRaw = haufePatternRaw(:);
result.haufePatternZ = haufePatternZ(:);
result.modelRelianceTAFC = modelRelianceTAFC;
result.featureMedian = prep.featureMedian(:)';
result.featureMean = prep.featureMean(:)';
result.featureScale = prep.featureScale(:)';
result.innerFoldIDByTrainSubject = foldID(:);
result.innerValidationLossMean = tune.meanLoss(:)';
result.innerValidationLossSE = tune.seLoss(:)';
result.lambdaPath = tune.lambdaPath(:)';

if numel(result.betaRaw) ~= nROI || numel(result.haufePatternRaw) ~= nROI
    error('Outer repeat %d produced unexpected coefficient or pattern sizes.', outerSplit.repeat);
end
end
