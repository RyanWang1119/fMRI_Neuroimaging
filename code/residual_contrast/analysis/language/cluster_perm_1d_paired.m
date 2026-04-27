function out = cluster_perm_1d_paired(X, Y, nPerm, alphaCluster)
% X, Y: [subjects x time]

if nargin < 3 || isempty(nPerm), nPerm = 5000; end
if nargin < 4 || isempty(alphaCluster), alphaCluster = 0.05; end

D = X - Y;
keep = all(isfinite(D), 2);
D = D(keep,:);

n = size(D,1);
T = size(D,2);

if n < 3
    error('Not enough subjects after removing NaNs.');
end

meanDiff = mean(D, 1);
sdDiff   = std(D, 0, 1);
tObs     = sqrt(n) * meanDiff ./ sdDiff;

tThresh = tinv(1 - alphaCluster/2, n - 1);

[posMassObs, posRanges] = get_cluster_masses(tObs,  tThresh,  1);
[negMassObs, negRanges] = get_cluster_masses(tObs, -tThresh, -1);

maxPosNull = zeros(nPerm,1);
maxNegNull = zeros(nPerm,1);

for p = 1:nPerm
    flips = 2*(rand(n,1) > 0.5) - 1;
    Dp = D .* flips;

    mp = mean(Dp,1);
    sp = std(Dp,0,1);
    tp = sqrt(n) * mp ./ sp;

    posMass = get_cluster_masses(tp,  tThresh,  1);
    negMass = get_cluster_masses(tp, -tThresh, -1);

    if ~isempty(posMass), maxPosNull(p) = max(posMass); end
    if ~isempty(negMass), maxNegNull(p) = max(negMass); end
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
        masses(i) = sum(tvals(idx));
    else
        masses(i) = sum(-tvals(idx));
    end
    ranges(i,:) = [starts(i), stops(i)];
end

masses = masses(:);
end