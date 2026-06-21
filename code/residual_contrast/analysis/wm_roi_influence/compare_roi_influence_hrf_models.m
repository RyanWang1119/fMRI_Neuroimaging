function comparisonTable = compare_roi_influence_hrf_models(resultFiles, outCsv)
%COMPARE_ROI_INFLUENCE_HRF_MODELS Descriptive parcel comparison across HRF models.
%
% This utility computes descriptive deltas only. It does not perform
% inferential tests for cross-model parcel comparisons because no paired
% resampling or permutation test is implemented here.
%
% resultFiles may be a struct with fields cHRF, cHRFderiv, and sHRF, or a
% 1x3 cell array in that order.

modelNames = {'cHRF', 'cHRFderiv', 'sHRF'};
files = resolve_files(resultFiles, modelNames);

D = cell(1, 3);
for m = 1:3
    S = load(files{m}, 'summaryTable');
    if ~isfield(S, 'summaryTable')
        error('File missing summaryTable: %s', files{m});
    end
    D{m} = S.summaryTable;
end

ref = D{1};
for m = 2:3
    if height(D{m}) ~= height(ref) || ~isequal(D{m}.roi_index, ref.roi_index)
        error('ROI rows do not match between %s and %s.', files{1}, files{m});
    end
end

comparisonTable = table();
comparisonTable.roi_index = ref.roi_index;
comparisonTable.roi_label = ref.roi_label;
comparisonTable.roi_network = ref.roi_network;

for m = 1:3
    suffix = modelNames{m};
    comparisonTable.(['meanAbsHaufe_' suffix]) = abs(D{m}.mean_haufe_raw);
    comparisonTable.(['selectionFrequency_' suffix]) = D{m}.selection_frequency;
    comparisonTable.(['modelReliance_' suffix]) = D{m}.mean_model_reliance_tafc;
end

comparisonTable.deltaMeanAbsHaufe_cHRFderiv_minus_cHRF = ...
    comparisonTable.meanAbsHaufe_cHRFderiv - comparisonTable.meanAbsHaufe_cHRF;
comparisonTable.deltaMeanAbsHaufe_sHRF_minus_cHRF = ...
    comparisonTable.meanAbsHaufe_sHRF - comparisonTable.meanAbsHaufe_cHRF;

comparisonTable.deltaSelectionFrequency_cHRFderiv_minus_cHRF = ...
    comparisonTable.selectionFrequency_cHRFderiv - comparisonTable.selectionFrequency_cHRF;
comparisonTable.deltaSelectionFrequency_sHRF_minus_cHRF = ...
    comparisonTable.selectionFrequency_sHRF - comparisonTable.selectionFrequency_cHRF;

comparisonTable.deltaModelReliance_cHRFderiv_minus_cHRF = ...
    comparisonTable.modelReliance_cHRFderiv - comparisonTable.modelReliance_cHRF;
comparisonTable.deltaModelReliance_sHRF_minus_cHRF = ...
    comparisonTable.modelReliance_sHRF - comparisonTable.modelReliance_cHRF;

if nargin >= 2 && ~isempty(outCsv)
    writetable(comparisonTable, outCsv);
end
end

function files = resolve_files(resultFiles, modelNames)
if isstruct(resultFiles)
    files = cell(1, 3);
    for m = 1:3
        if ~isfield(resultFiles, modelNames{m})
            error('resultFiles missing field %s.', modelNames{m});
        end
        files{m} = resultFiles.(modelNames{m});
    end
elseif iscell(resultFiles) && numel(resultFiles) == 3
    files = resultFiles;
else
    error('resultFiles must be a struct with cHRF/cHRFderiv/sHRF fields or a 1x3 cell array.');
end

for m = 1:3
    if ~isfile(files{m})
        error('Result file not found: %s', files{m});
    end
end
end
