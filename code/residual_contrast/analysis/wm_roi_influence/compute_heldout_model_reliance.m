function [deltaTAFC, tafcAfterSwaps] = compute_heldout_model_reliance(DtestZ, betaStandardized, selectedMask, cfg, seed)
%COMPUTE_HELDOUT_MODEL_RELIANCE Paired half-swap held-out reliance.
%
% This is held-out model reliance, not a causal measure and not total
% biological information in a parcel. Correlated parcels can substitute
% for each other, so a small performance drop does not imply absence of
% task-related information.

DtestZ = double(DtestZ);
betaStandardized = double(betaStandardized(:));
selectedMask = logical(selectedMask(:));

[nTest, nROI] = size(DtestZ);
if numel(betaStandardized) ~= nROI
    error('betaStandardized length must match DtestZ columns.');
end
if numel(selectedMask) ~= nROI
    error('selectedMask length must match DtestZ columns.');
end

nSwaps = cfg.importance.nRandomSwaps;
deltaTAFC = zeros(nROI, nSwaps, 'single');
tafcAfterSwaps = nan(nROI, nSwaps, 'single');

baselineMargins = DtestZ * betaStandardized;
baselineTAFC = mean(baselineMargins > 0, 'omitnan');

if cfg.importance.selectedOnly
    parcelsToTest = find(selectedMask);
else
    parcelsToTest = (1:nROI)';
end

if isempty(parcelsToTest) || nTest < 1
    return
end

stream = RandStream('mt19937ar', 'Seed', seed);
nSwapSubjects = max(1, floor(nTest / 2));

for p = parcelsToTest(:)'
    for s = 1:nSwaps
        swapIdx = randperm(stream, nTest, nSwapSubjects);
        Dswap = DtestZ;
        Dswap(swapIdx, p) = -Dswap(swapIdx, p);
        tafc = mean((Dswap * betaStandardized) > 0, 'omitnan');
        tafcAfterSwaps(p, s) = single(tafc);
        deltaTAFC(p, s) = single(baselineTAFC - tafc);
    end
end
end
