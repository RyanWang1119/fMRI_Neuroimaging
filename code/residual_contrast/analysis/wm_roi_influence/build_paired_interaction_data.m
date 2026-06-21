function pairedData = build_paired_interaction_data(residCurves, cfg, windowTRs, conditionMap)
%BUILD_PAIRED_INTERACTION_DATA Compute A, B, and D subject-level WM contrasts.
%
% For each subject and selected post-stimulus window:
%   A_i = mean_t residual(2-back Category A) - residual(0-back Category A)
%   B_i = mean_t residual(2-back Category B) - residual(0-back Category B)
%   D_i = A_i - B_i
%
% A, B, and D are [nSubject x nROI]. Downstream logistic models use the
% symmetric paired representation [D/2; -D/2] and never independently
% shuffle task-condition rows.

if ndims(residCurves) ~= 4
    error('residCurves must be [nROI x nTime x nSubject x nCondition].');
end

[nROI, nTime, nSubject, nCondition] = size(residCurves);
if nCondition ~= numel(conditionMap.conditionLabels)
    error('Expected %d conditions, found %d.', numel(conditionMap.conditionLabels), nCondition);
end
if any(windowTRs < 1) || any(windowTRs > nTime)
    error('windowTRs contains indices outside 1..%d.', nTime);
end

contrastIdx = resolve_contrast_index(cfg, conditionMap);
catPair = conditionMap.interactionCategoryPairs(contrastIdx, :);

catA = catPair(1);
catB = catPair(2);

idx0A = conditionMap.zeroBackIdx(catA);
idx2A = conditionMap.twoBackIdx(catA);
idx0B = conditionMap.zeroBackIdx(catB);
idx2B = conditionMap.twoBackIdx(catB);

A = mean_window_difference(residCurves, windowTRs, idx2A, idx0A);
B = mean_window_difference(residCurves, windowTRs, idx2B, idx0B);
D = A - B;

pairedData = struct();
pairedData.A = A;
pairedData.B = B;
pairedData.D = D;
pairedData.contrastIndex = contrastIdx;
pairedData.contrastLabel = conditionMap.interactionLabels{contrastIdx};
pairedData.categoryA = conditionMap.categoryNames{catA};
pairedData.categoryB = conditionMap.categoryNames{catB};
pairedData.conditionIndices = [idx2A idx0A idx2B idx0B];
pairedData.windowTRs = windowTRs(:)';
pairedData.nROI = nROI;
pairedData.nSubject = nSubject;
end

function X = mean_window_difference(R, windowTRs, idx2, idx0)
diffCurve = R(:, windowTRs, :, idx2) - R(:, windowTRs, :, idx0);
XroiSubj = squeeze(mean(diffCurve, 2, 'omitnan')); % [nROI x nSubject]
if isvector(XroiSubj)
    XroiSubj = reshape(XroiSubj, size(R,1), []);
end
X = XroiSubj'; % [nSubject x nROI]
end

function contrastIdx = resolve_contrast_index(cfg, conditionMap)
if isfield(cfg, 'contrastIndex') && ~isempty(cfg.contrastIndex)
    contrastIdx = cfg.contrastIndex;
    if ~isscalar(contrastIdx) || contrastIdx < 1 || contrastIdx > numel(conditionMap.interactionLabels)
        error('cfg.contrastIndex must be an integer from 1 to %d.', numel(conditionMap.interactionLabels));
    end
    contrastIdx = round(contrastIdx);
    return
end

if ~isfield(cfg, 'contrastName') || isempty(cfg.contrastName)
    error('Specify cfg.contrastName or cfg.contrastIndex.');
end

q = normalize_label(cfg.contrastName);
labels = cellfun(@normalize_label, conditionMap.interactionLabels, 'UniformOutput', false);
contrastIdx = find(strcmp(q, labels), 1);
if ~isempty(contrastIdx)
    return
end

% Friendly aliases such as "body vs faces".
for i = 1:size(conditionMap.interactionCategoryPairs, 1)
    pair = conditionMap.interactionCategoryPairs(i, :);
    aAliases = conditionMap.categoryAliases{pair(1)};
    bAliases = conditionMap.categoryAliases{pair(2)};
    for ia = 1:numel(aAliases)
        for ib = 1:numel(bAliases)
            aliases = { ...
                sprintf('%s_vs_%s', aAliases{ia}, bAliases{ib}), ...
                sprintf('%s vs %s', aAliases{ia}, bAliases{ib}), ...
                sprintf('%s-loaddiff-vs-%s-loaddiff', aAliases{ia}, bAliases{ib})};
            aliases = cellfun(@normalize_label, aliases, 'UniformOutput', false);
            if any(strcmp(q, aliases))
                contrastIdx = i;
                return
            end
        end
    end
end

error('Could not resolve contrastName "%s". Available contrasts: %s', ...
    cfg.contrastName, strjoin(conditionMap.interactionLabels(:)', ', '));
end

function y = normalize_label(x)
y = lower(char(x));
y = strrep(y, 'faces', 'face');
y = strrep(y, 'places', 'place');
y = strrep(y, 'tools', 'tool');
y = regexprep(y, '[^a-z0-9]+', '_');
y = regexprep(y, '_+', '_');
y = regexprep(y, '^_|_$', '');
end
