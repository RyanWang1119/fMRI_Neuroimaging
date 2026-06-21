function signFlipResults = run_paired_signflip_maxT(D, cfg)
%RUN_PAIRED_SIGNFLIP_MAXT Population-level paired sign-flip maxT test.
%
% This non-machine-learning regional test uses the unstandardized
% subject-level window-averaged paired difference D_i = A_i - B_i. It does
% not use decoder coefficients, feature selections, or held-out performance.

D = double(D);
[nSubject, nROI] = size(D);
nPerm = cfg.permutation.nPerm;
alpha = cfg.permutation.alpha;

[tObs, meanRawDifference, effectSizeCohenDz, nValid] = compute_t_by_parcel(D);
maxAbsT_perm = nan(nPerm, 1);

stream = RandStream('mt19937ar', 'Seed', cfg.seed + 7919);
for p = 1:nPerm
    signs = 2 * (rand(stream, nSubject, 1) > 0.5) - 1;
    Dperm = bsxfun(@times, D, signs);
    tPerm = compute_t_by_parcel(Dperm);
    maxAbsT_perm(p) = max(abs(tPerm), [], 'omitnan');
end

pFWER = nan(nROI, 1);
for r = 1:nROI
    if isfinite(tObs(r))
        pFWER(r) = (1 + sum(maxAbsT_perm >= abs(tObs(r)))) / (nPerm + 1);
    end
end
significantFWER = pFWER < alpha;

signFlipResults = struct();
signFlipResults.tObs = tObs(:);
signFlipResults.pFWER = pFWER(:);
signFlipResults.significantFWER = significantFWER(:);
signFlipResults.meanRawDifference = meanRawDifference(:);
signFlipResults.effectSizeCohenDz = effectSizeCohenDz(:);
signFlipResults.nValid = nValid(:);
signFlipResults.maxAbsT_perm = maxAbsT_perm(:);
signFlipResults.nPerm = nPerm;
signFlipResults.alpha = alpha;
signFlipResults.seed = cfg.seed + 7919;
signFlipResults.note = ['Subject-level paired sign flips with max absolute t-statistic ', ...
    'across parcels; independent of decoder outputs.'];
end

function [tStat, meanD, dz, nValid] = compute_t_by_parcel(D)
nROI = size(D, 2);
tStat = nan(nROI, 1);
meanD = nan(nROI, 1);
dz = nan(nROI, 1);
nValid = zeros(nROI, 1);

for r = 1:nROI
    x = D(:, r);
    x = x(isfinite(x));
    n = numel(x);
    nValid(r) = n;
    if n < 2
        continue
    end
    m = mean(x);
    s = std(x, 0);
    meanD(r) = m;
    if s <= eps
        if abs(m) <= eps
            tStat(r) = 0;
            dz(r) = 0;
        else
            tStat(r) = sign(m) * Inf;
            dz(r) = sign(m) * Inf;
        end
    else
        tStat(r) = m / (s / sqrt(n));
        dz(r) = m / s;
    end
end
end
