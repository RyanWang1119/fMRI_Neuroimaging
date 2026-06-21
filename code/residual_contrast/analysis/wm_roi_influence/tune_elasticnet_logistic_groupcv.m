function tune = tune_elasticnet_logistic_groupcv(Dtrain, cfg, foldID)
%TUNE_ELASTICNET_LOGISTIC_GROUPCV Custom grouped CV for lassoglm.
%
% The lambda path is constructed from the full outer-training set only.
% Each inner fold refits preprocessing on inner-training subjects and
% evaluates negative logistic log-likelihood on held-out subjects, keeping
% both synthetic paired rows from a subject in the same fold.

if exist('lassoglm', 'file') ~= 2
    error('lassoglm is required for elastic-net logistic regression.');
end

Dtrain = double(Dtrain);
nTrain = size(Dtrain, 1);
if numel(foldID) ~= nTrain
    error('foldID length must match number of training subjects.');
end
foldID = foldID(:);
folds = unique(foldID(:))';
if numel(folds) < 2
    error('Need at least two inner folds.');
end

outerPrep = apply_training_preprocessing(Dtrain, []);
ypairOuter = paired_labels(nTrain);
[~, fitInfoPath] = lassoglm(outerPrep.XpairTrainZ, ypairOuter, 'binomial', ...
    'Alpha', cfg.model.alpha, ...
    'NumLambda', cfg.model.numLambda, ...
    'Standardize', false);
lambdaPath = fitInfoPath.Lambda(:)';
lambdaPath = sort(unique(lambdaPath), 'descend');
if isempty(lambdaPath)
    error('lassoglm returned an empty lambda path.');
end

loss = nan(numel(folds), numel(lambdaPath));

for f = 1:numel(folds)
    valMask = (foldID == folds(f));
    trainMask = ~valMask;
    if ~any(trainMask) || ~any(valMask)
        error('Inner fold %d has empty train or validation set.', folds(f));
    end

    prep = apply_training_preprocessing(Dtrain(trainMask, :), Dtrain(valMask, :));
    yInner = paired_labels(sum(trainMask));
    yVal = paired_labels(sum(valMask));

    [B, fitInfo] = lassoglm(prep.XpairTrainZ, yInner, 'binomial', ...
        'Alpha', cfg.model.alpha, ...
        'Lambda', lambdaPath, ...
        'Standardize', false);

    for l = 1:numel(fitInfo.Lambda)
        lambdaIdx = find(abs(lambdaPath - fitInfo.Lambda(l)) <= max(eps(lambdaPath), eps(fitInfo.Lambda(l))), 1);
        if isempty(lambdaIdx)
            continue
        end
        eta = prep.XpairApplyZ * B(:, l) + fitInfo.Intercept(l);
        loss(f, lambdaIdx) = logistic_nll(yVal, eta);
    end
end

meanLoss = mean(loss, 1, 'omitnan');
nValid = sum(isfinite(loss), 1);
seLoss = std(loss, 0, 1, 'omitnan') ./ sqrt(max(nValid, 1));

if all(~isfinite(meanLoss))
    error('All inner validation losses are non-finite.');
end

[minLoss, minIdx] = min(meanLoss);
if cfg.model.useOneSE
    threshold = minLoss + seLoss(minIdx);
    candidateIdx = find(meanLoss <= threshold & isfinite(meanLoss));
    [~, local] = max(lambdaPath(candidateIdx));
    selectedIdx = candidateIdx(local);
else
    selectedIdx = minIdx;
end

tune = struct();
tune.lambdaPath = lambdaPath;
tune.loss = loss;
tune.meanLoss = meanLoss;
tune.seLoss = seLoss;
tune.minLoss = minLoss;
tune.minLambda = lambdaPath(minIdx);
tune.selectedLambda = lambdaPath(selectedIdx);
tune.selectedLambdaIndex = selectedIdx;
tune.foldID = foldID;
tune.rule = ternary(cfg.model.useOneSE, 'one-standard-error largest lambda', 'minimum mean validation loss');
end

function y = paired_labels(nSubject)
y = [ones(nSubject, 1); zeros(nSubject, 1)];
end

function val = logistic_nll(y, eta)
eta = double(eta(:));
y = double(y(:));
val = mean(max(eta, 0) - y .* eta + log1p(exp(-abs(eta))), 'omitnan');
end

function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end
