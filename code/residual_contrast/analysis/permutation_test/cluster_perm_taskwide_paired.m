function Results = cluster_perm_taskwide_paired(Xcell, Ycell, nPerm, alphaCluster)
% cluster_perm_taskwide_paired
% Task-wide within-subject cluster permutation test across time and contrasts.
%
% INPUT
%   Xcell, Ycell : cell arrays, one cell per contrast
%                  each cell is [subjects x time]
%   nPerm        : number of sign-flip permutations
%   alphaCluster : cluster-forming alpha, typically 0.05
%
% OUTPUT
%   Results(k) for each contrast k:
%       .meanDiff
%       .tObs
%       .nSubjects
%       .posClusters
%       .negClusters

if nargin < 3 || isempty(nPerm)
    nPerm = 5000;
end
if nargin < 4 || isempty(alphaCluster)
    alphaCluster = 0.05;
end

K = numel(Xcell);
assert(K == numel(Ycell), 'Xcell and Ycell must have same number of contrasts.');

S = size(Xcell{1},1);
T = size(Xcell{1},2);

Dfull   = cell(K,1);
keepRow = cell(K,1);
nVec    = zeros(K,1);
tThresh = zeros(K,1);

meanDiffCell = cell(K,1);
tObsCell     = cell(K,1);
posMassObs   = cell(K,1);
posRangesObs = cell(K,1);
negMassObs   = cell(K,1);
negRangesObs = cell(K,1);

for k = 1:K
    X = Xcell{k};
    Y = Ycell{k};

    assert(all(size(X) == size(Y)), 'Contrast %d: X and Y size mismatch.', k);
    assert(size(X,1) == S && size(X,2) == T, 'Contrast %d: inconsistent matrix size.', k);

    D = X - Y;                     % [S x T]
    keep = all(isfinite(D), 2);    % complete-case subjects for this contrast
    Dk = D(keep, :);

    n = size(Dk,1);
    if n < 3
        error('Contrast %d has fewer than 3 usable subjects.', k);
    end

    md = mean(Dk, 1);
    sd = std(Dk, 0, 1);

    tObs = zeros(1, T);
    ok = sd > 0;
    tObs(ok) = sqrt(n) * md(ok) ./ sd(ok);

    thr = tinv(1 - alphaCluster/2, n - 1);

    [pm, pr] = get_cluster_masses(tObs,  thr,  1);
    [nm, nr] = get_cluster_masses(tObs, -thr, -1);

    Dfull{k}      = D;
    keepRow{k}    = keep;
    nVec(k)       = n;
    tThresh(k)    = thr;
    meanDiffCell{k} = md;
    tObsCell{k}     = tObs;
    posMassObs{k}   = pm;
    posRangesObs{k} = pr;
    negMassObs{k}   = nm;
    negRangesObs{k} = nr;
end

% Task-wide null: max cluster mass across ALL contrasts in the task
maxPosNull = zeros(nPerm,1);
maxNegNull = zeros(nPerm,1);

for p = 1:nPerm
    signFull = 2 * (rand(S,1) > 0.5) - 1;   % one sign per subject, shared across contrasts

    curMaxPos = 0;
    curMaxNeg = 0;

    for k = 1:K
        Dp = Dfull{k}(keepRow{k}, :) .* signFull(keepRow{k});

        mp = mean(Dp, 1);
        sp = std(Dp, 0, 1);

        tp = zeros(1, T);
        ok = sp > 0;
        tp(ok) = sqrt(nVec(k)) * mp(ok) ./ sp(ok);

        posMass = get_cluster_masses(tp,  tThresh(k),  1);
        negMass = get_cluster_masses(tp, -tThresh(k), -1);

        if ~isempty(posMass)
            curMaxPos = max(curMaxPos, max(posMass));
        end
        if ~isempty(negMass)
            curMaxNeg = max(curMaxNeg, max(negMass));
        end
    end

    maxPosNull(p) = curMaxPos;
    maxNegNull(p) = curMaxNeg;
end

Results = struct([]);
for k = 1:K
    Results(k).meanDiff  = meanDiffCell{k};
    Results(k).tObs      = tObsCell{k};
    Results(k).nSubjects = nVec(k);

    Results(k).posClusters = struct('start',{},'stop',{},'mass',{},'p',{});
    for j = 1:numel(posMassObs{k})
        pval = (1 + sum(maxPosNull >= posMassObs{k}(j))) / (nPerm + 1);
        Results(k).posClusters(j).start = posRangesObs{k}(j,1);
        Results(k).posClusters(j).stop  = posRangesObs{k}(j,2);
        Results(k).posClusters(j).mass  = posMassObs{k}(j);
        Results(k).posClusters(j).p     = pval;
    end

    Results(k).negClusters = struct('start',{},'stop',{},'mass',{},'p',{});
    for j = 1:numel(negMassObs{k})
        pval = (1 + sum(maxNegNull >= negMassObs{k}(j))) / (nPerm + 1);
        Results(k).negClusters(j).start = negRangesObs{k}(j,1);
        Results(k).negClusters(j).stop  = negRangesObs{k}(j,2);
        Results(k).negClusters(j).mass  = negMassObs{k}(j);
        Results(k).negClusters(j).p     = pval;
    end
end

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
        masses(i,1) = sum(tvals(idx));    %#ok<AGROW>
    else
        masses(i,1) = sum(-tvals(idx));   %#ok<AGROW>
    end
    ranges(i,:) = [starts(i), stops(i)];  %#ok<AGROW>
end
end