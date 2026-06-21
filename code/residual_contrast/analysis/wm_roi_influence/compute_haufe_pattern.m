function [haufePatternZ, haufePatternRaw] = compute_haufe_pattern(XpairTrainZ, betaStandardized, featureScale)
%COMPUTE_HAUFE_PATTERN Haufe pattern for standardized paired design.
%
% Haufe pattern:
%   SigmaX = cov(Xpair_train_z);
%   haufePatternZ = SigmaX * betaStandardized;
%
% Raw-unit convention:
%   featureScale is the training scale used to standardize the synthetic
%   paired rows [D/2; -D/2]. Multiplying haufePatternZ by featureScale
%   expresses the pattern in raw paired-row residual-difference units.

betaStandardized = betaStandardized(:);
featureScale = featureScale(:);

if size(XpairTrainZ, 2) ~= numel(betaStandardized)
    error('XpairTrainZ columns must match betaStandardized length.');
end

SigmaX = cov(XpairTrainZ);
haufePatternZ = SigmaX * betaStandardized;
haufePatternRaw = featureScale .* haufePatternZ;
end
