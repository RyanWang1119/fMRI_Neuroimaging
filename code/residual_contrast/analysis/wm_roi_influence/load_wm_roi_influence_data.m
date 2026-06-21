function [residCurves, subjectIDs, roiLabels, metadata] = load_wm_roi_influence_data(cfg)
%LOAD_WM_ROI_INFLUENCE_DATA Load WM HRF-model residual curves.
%
% Output residCurves is [nROI x nTime x nSubject x nCondition].

if ~isfield(cfg, 'hrfModelName') || isempty(cfg.hrfModelName)
    error('cfg.hrfModelName is required.');
end

dataFile = resolve_hrf_model_file(cfg);
vars = whos('-file', dataFile);
varNames = {vars.name};

if ismember('Results', varNames)
    dataVariable = 'Results';
else
    candidates = varNames([vars.bytes] > 0);
    dataVariable = '';
    for i = 1:numel(candidates)
        vi = vars(strcmp(varNames, candidates{i}));
        if numel(vi.size) >= 4 && vi.size(end) == 8
            dataVariable = candidates{i};
            break
        end
    end
    if isempty(dataVariable)
        error('Could not find WM residual array variable in %s. Expected variable Results.', dataFile);
    end
end

S = load(dataFile, dataVariable);
residCurves = S.(dataVariable);
residCurves = canonicalize_wm_array(residCurves, dataFile);

[nROI, nTime, nSubject, nCondition] = size(residCurves);
if nROI ~= 489 || nTime ~= 41 || nCondition ~= 8
    error('Unexpected WM residual size in %s: [%s]. Expected [489 x 41 x nSubject x 8].', ...
        dataFile, num2str(size(residCurves)));
end

[subjectIDs, subjectIDSource] = get_subject_ids(cfg, dataFile, varNames, nSubject);
[roiLabels, roiLabelSource] = get_roi_labels(cfg, dataFile, varNames, nROI);

metadata = struct();
metadata.dataFile = dataFile;
metadata.dataVariable = dataVariable;
metadata.dataSize = [nROI nTime nSubject nCondition];
metadata.subjectIDSource = subjectIDSource;
metadata.roiLabelSource = roiLabelSource;
metadata.importanceComputation = ['Held-out model reliance is computed by paired half-swaps ', ...
    'for selected parcels only by default; nonselected parcels receive zero deltaTAFC.'];
end

function dataFile = resolve_hrf_model_file(cfg)
if isfield(cfg, 'dataFile') && ~isempty(cfg.dataFile)
    dataFile = cfg.dataFile;
else
    if ~isfield(cfg, 'dataDir') || isempty(cfg.dataDir)
        if isfield(cfg, 'repoRoot') && ~isempty(cfg.repoRoot)
            cfg.dataDir = fullfile(cfg.repoRoot, 'data', 'task_residual');
        else
            cfg.dataDir = fullfile(pwd, 'data', 'task_residual');
        end
    end
    switch lower(char(cfg.hrfModelName))
        case 'chrf'
            fname = 'WMcHRF.mat';
        case 'chrfderiv'
            fname = 'WMcHRFderiv.mat';
        case 'shrf'
            fname = 'WMsHRF.mat';
        otherwise
            error('Unsupported cfg.hrfModelName: %s. Use cHRF, cHRFderiv, or sHRF.', cfg.hrfModelName);
    end
    dataFile = fullfile(cfg.dataDir, fname);
end

if ~isfile(dataFile)
    error('WM residual data file not found: %s', dataFile);
end
end

function R = canonicalize_wm_array(Rin, dataFile)
sz = size(Rin);
sz(end+1:4) = 1;

if sz(1) == 489 && sz(2) == 41 && sz(4) == 8
    R = Rin;
elseif sz(2) == 489 && sz(3) == 41 && sz(4) == 8
    R = permute(Rin, [2 3 1 4]);
else
    error('Unexpected WM array format in %s: [%s]. Expected [489 x 41 x S x 8] or [S x 489 x 41 x 8].', ...
        dataFile, num2str(sz));
end
end

function [subjectIDs, source] = get_subject_ids(cfg, dataFile, varNames, nSubject)
if isfield(cfg, 'subjectIDs') && ~isempty(cfg.subjectIDs)
    subjectIDs = cfg.subjectIDs(:);
    source = 'cfg.subjectIDs';
elseif ismember('subjectIDs', varNames)
    X = load(dataFile, 'subjectIDs');
    subjectIDs = X.subjectIDs(:);
    source = sprintf('%s:subjectIDs', dataFile);
elseif ismember('subject_ids', varNames)
    X = load(dataFile, 'subject_ids');
    subjectIDs = X.subject_ids(:);
    source = sprintf('%s:subject_ids', dataFile);
elseif ismember('subjIDs', varNames)
    X = load(dataFile, 'subjIDs');
    subjectIDs = X.subjIDs(:);
    source = sprintf('%s:subjIDs', dataFile);
else
    subjectIDs = (1:nSubject)';
    source = 'ordinal row index; HRF-model files do not store real HCP subject IDs';
end

if numel(subjectIDs) ~= nSubject
    error('subjectIDs has %d entries, but data has %d subjects.', numel(subjectIDs), nSubject);
end
end

function [roiLabels, source] = get_roi_labels(cfg, dataFile, varNames, nROI)
if isfield(cfg, 'roiLabels') && ~isempty(cfg.roiLabels)
    roiLabels = normalize_roi_labels(cfg.roiLabels, nROI);
    source = 'cfg.roiLabels';
    return
end

candidateNames = {'roiLabels', 'region_names', 'RegionNames', 'labels'};
for i = 1:numel(candidateNames)
    if ismember(candidateNames{i}, varNames)
        X = load(dataFile, candidateNames{i});
        roiLabels = normalize_roi_labels(X.(candidateNames{i}), nROI);
        source = sprintf('%s:%s', dataFile, candidateNames{i});
        return
    end
end

if exist('load_atlas', 'file') == 2
    try
        atlasObj = load_atlas('canlab2018');
        atlasLabels = [];
        if isstruct(atlasObj) && isfield(atlasObj, 'labels')
            atlasLabels = atlasObj.labels;
        elseif isobject(atlasObj) && isprop(atlasObj, 'labels')
            atlasLabels = atlasObj.labels;
        end
        if ~isempty(atlasLabels) && numel(atlasLabels) == nROI
            roiLabels = normalize_roi_labels(atlasLabels, nROI);
            source = 'load_atlas(''canlab2018'').labels';
            return
        end
    catch
        % Fall through to generic labels.
    end
end

roiLabels = normalize_roi_labels(arrayfun(@(r) sprintf('Parcel_%03d', r), ...
    (1:nROI)', 'UniformOutput', false), nROI);
source = 'generic Parcel_### labels; canlab2018 atlas labels unavailable';
end

function T = normalize_roi_labels(labelsIn, nROI)
if istable(labelsIn)
    T = labelsIn;
    if height(T) ~= nROI
        error('ROI label table has %d rows; expected %d.', height(T), nROI);
    end
    if ~ismember('roi_index', T.Properties.VariableNames)
        T.roi_index = (1:nROI)';
    end
    if ~ismember('roi_label', T.Properties.VariableNames)
        firstText = find(varfun(@(x) iscellstr(x) || isstring(x) || iscategorical(x), ...
            T, 'OutputFormat', 'uniform'), 1);
        if isempty(firstText)
            T.roi_label = arrayfun(@(r) sprintf('Parcel_%03d', r), (1:nROI)', 'UniformOutput', false);
        else
            T.roi_label = cellstr(T{:, firstText});
        end
    end
    if ~ismember('roi_network', T.Properties.VariableNames)
        T.roi_network = repmat({''}, nROI, 1);
    end
    T = T(:, {'roi_index', 'roi_label', 'roi_network'});
    return
end

if isstring(labelsIn) || iscategorical(labelsIn)
    labels = cellstr(labelsIn(:));
elseif iscell(labelsIn)
    labels = labelsIn(:);
else
    labels = arrayfun(@(r) sprintf('Parcel_%03d', r), (1:nROI)', 'UniformOutput', false);
end

if numel(labels) ~= nROI
    error('ROI labels have %d entries; expected %d.', numel(labels), nROI);
end

T = table((1:nROI)', labels(:), repmat({''}, nROI, 1), ...
    'VariableNames', {'roi_index', 'roi_label', 'roi_network'});
end
