function prep = apply_training_preprocessing(Dtrain, Dapply)
%APPLY_TRAINING_PREPROCESSING Fit train-only imputation and scaling.
%
% Dtrain and Dapply are subject-level paired differences [subject x ROI].
% Missing D values are imputed with training-feature medians only.
%
% The logistic design is pair-centered:
%   Xpair = [ D/2 ; -D/2 ].
% Means and scales are fit on the paired training rows. Because the paired
% rows are symmetric after imputation, the feature mean is zero apart from
% numerical precision. For held-out two-alternative margins, the paired-row
% centering cancels, so D_z is D_imputed ./ featureScale.

if nargin < 2
    Dapply = [];
end

Dtrain = double(Dtrain);
Dapply = double(Dapply);
if size(Dtrain, 1) < 1
    error('Dtrain must contain at least one subject.');
end

nFeature = size(Dtrain, 2);
if ~isempty(Dapply) && size(Dapply, 2) ~= nFeature
    error('Dapply has %d features; expected %d.', size(Dapply, 2), nFeature);
end

featureMedian = median(Dtrain, 1, 'omitnan');
allMissing = isnan(featureMedian);
featureMedian(allMissing) = 0;

DtrainImp = impute_with_values(Dtrain, featureMedian);
DapplyImp = impute_with_values(Dapply, featureMedian);

XpairTrainRaw = [DtrainImp / 2; -DtrainImp / 2];
featureMean = mean(XpairTrainRaw, 1, 'omitnan');
Xcenter = XpairTrainRaw - featureMean;
featureScale = std(Xcenter, 0, 1, 'omitnan');
badScale = ~isfinite(featureScale) | featureScale <= eps;
featureScale(badScale) = 1;

XpairTrainZ = Xcenter ./ featureScale;
if isempty(DapplyImp)
    XpairApplyZ = [];
    DapplyZ = [];
else
    XpairApplyRaw = [DapplyImp / 2; -DapplyImp / 2];
    XpairApplyZ = (XpairApplyRaw - featureMean) ./ featureScale;
    DapplyZ = DapplyImp ./ featureScale;
end

prep = struct();
prep.featureMedian = featureMedian;
prep.featureMean = featureMean;
prep.featureScale = featureScale;
prep.allMissingFeature = allMissing;
prep.constantFeature = badScale;
prep.DtrainImputed = DtrainImp;
prep.DapplyImputed = DapplyImp;
prep.XpairTrainZ = XpairTrainZ;
prep.XpairApplyZ = XpairApplyZ;
prep.DtrainZ = DtrainImp ./ featureScale;
prep.DapplyZ = DapplyZ;
end

function X = impute_with_values(X, values)
if isempty(X)
    return
end
nanMask = isnan(X);
if ~any(nanMask(:))
    return
end
for j = 1:size(X, 2)
    bad = nanMask(:, j);
    if any(bad)
        X(bad, j) = values(j);
    end
end
end
