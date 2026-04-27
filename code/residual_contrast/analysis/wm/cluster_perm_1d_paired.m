function out = cluster_perm_1d_paired(X, Y, nPerm, alphaCluster)
% cluster_perm_1d_paired
% Within-subject sign-flip cluster permutation test across time.
%
% INPUT
%   X, Y         : [subjects x time]
%   nPerm        : number of sign-flip permutations
%   alphaCluster : cluster-forming alpha, typically 0.05
%
% OUTPUT
%   out.meanDiff
%   out.tObs
%   out.posClusters
%   out.negClusters

if nargin < 3 || isempty(nPerm)
    nPerm = 5000;
end
if nargin < 4 || isempty(alphaCluster)
    alphaCluster = 0.05;
end

assert(ismatrix(X) && ismatrix(Y), 'X and Y must be 2D [subjects x time].');
assert(all(size(X) == size(Y)), 'X and Y must have the same size.');

D = X - Y;                         % paired difference curves
keep = all(isfinite(D), 2);        % keep subjects with finite values at all time points
D = D(keep,:);

n = size(D,1);
T = size(D,2);

if n < 3
    error('Not enough subjects after removing NaN rows.');
end

meanDiff = mean(D, 1);
sdDiff   = std(D, 0, 1);

tObs = nan(1, T);
ok = sdDiff > 0;
tObs(ok) = sqrt(n) * meanDiff(ok) ./ sdDiff(ok);
tObs(~ok) = 0;

tThresh = tinv(1 - alphaCluster/2, n - 1);

[posMassObs, posRanges] = get_cluster_masses(tObs,  tThresh,  1);
[negMassObs, negRanges] = get_cluster_masses(tObs, -tThresh, -1);

maxPosNull = zeros(nPerm,1);
maxNegNull = zeros(nPerm,1);

for p = 1:nPerm
    flips = 2 * (rand(n,1) > 0.5) - 1;    % +/- 1
    Dp = D .* flips;

    mp = mean(Dp, 1);
    sp = std(Dp, 0, 1);

    tp = zeros(1, T);
    okp = sp > 0;
    tp(okp) = sqrt(n) * mp(okp) ./ sp(okp);

    posMass = get_cluster_masses(tp,  tThresh,  1);
    negMass = get_cluster_masses(tp, -tThresh, -1);

    if ~isempty(posMass)
        maxPosNull(p) = max(posMass);
    end
    if ~isempty(negMass)
        maxNegNull(p) = max(negMass);
    end
end

posClusters = struct('start',{},'stop',{},'mass',{},'p',{});
for j = 1:numel(posMassObs)
    pval = (1 + sum(maxPosNull >= posMassObs(j))) / (nPerm + 1);
    posClusters(j).start = posRanges(j,1);
    posClusters(j).stop  = posRanges(j,2);
    posClusters(j).mass  = posMassObs(j);
    posClusters(j).p     = pval;
end

negClusters = struct('start',{},'stop',{},'mass',{},'p',{});
for j = 1:numel(negMassObs)
    pval = (1 + sum(maxNegNull >= negMassObs(j))) / (nPerm + 1);
    negClusters(j).start = negRanges(j,1);
    negClusters(j).stop  = negRanges(j,2);
    negClusters(j).mass  = negMassObs(j);
    negClusters(j).p     = pval;
end

out.meanDiff    = meanDiff;
out.tObs        = tObs;
out.posClusters = posClusters;
out.negClusters = negClusters;
end

function [masses, ranges] = get_cluster_masses(tvals, thresh, modeFlag)
% modeFlag = 1 for positive clusters, -1 for negative clusters

if modeFlag == 1
    mask = tvals > thresh;
else
    mask = tvals < thresh;
end

d = diff([0 mask 0]);
starts = find(d == 1);
stops  = find(d == -1) - 1;

masses = [];
ranges = [];

for i = 1:numel(starts)
    idx = starts(i):stops(i);
    if modeFlag == 1
        masses(i,1) = sum(tvals(idx));   %#ok<AGROW>
    else
        masses(i,1) = sum(-tvals(idx));  %#ok<AGROW>
    end
    ranges(i,:) = [starts(i), stops(i)]; %#ok<AGROW>
end
end