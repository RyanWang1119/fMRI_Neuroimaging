function out = run_roi_influence_elasticnet(cfg)
% RUN_ROI_INFLUENCE_ELASTICNET
%
% Separate WM residual-decoding module for stable, high-confidence
% influential candidate parcels using nested elastic-net logistic regression.
%
% Repository conventions discovered from the existing WM residual-decoding
% code before this module was written:
%   - HRF-model files data/task_residual/WMcHRF.mat, WMcHRFderiv.mat, and
%     WMsHRF.mat contain variable Results with size [489 x 41 x 410 x 8].
%   - The older raw file HRF_Resid_WMblock_LR.mat contains Results with size
%     [489 x 41 x 393 x 8] plus Flags, and the current SVM wrapper
%     canonicalizes the same residual-FIR order to [subjects x ROI x time x condition].
%   - WM condition slots are determined from functional_clustering/build_wm_contrast.m
%     and code/residual_contrast/data_extract/PredictContrast.m:
%       1=0bk_body,  2=0bk_faces,  3=0bk_places,  4=0bk_tools,
%       5=2bk_body,  6=2bk_faces,  7=2bk_places,  8=2bk_tools.
%   - The six working-memory interaction contrasts are:
%       [5 1 6 2] Body_LoadDiff_vs_Face_LoadDiff,
%       [5 1 7 3] Body_LoadDiff_vs_Place_LoadDiff,
%       [5 1 8 4] Body_LoadDiff_vs_Tool_LoadDiff,
%       [6 2 7 3] Face_LoadDiff_vs_Place_LoadDiff,
%       [6 2 8 4] Face_LoadDiff_vs_Tool_LoadDiff,
%       [7 3 8 4] Place_LoadDiff_vs_Tool_LoadDiff.
%   - Existing WM plotting code uses TR = 0.72 seconds.
%   - Real HCP subject IDs are not stored in the three HRF-model .mat files;
%     this module uses ordinal subject IDs unless cfg.subjectIDs is supplied.
%
% cfg.windowSec is required. It must be pre-specified from prior reasoning
% or selected using training data only; this driver does not pick a window
% from the observed full-dataset peak.

if nargin < 1
    cfg = struct();
end
if ~isstruct(cfg)
    error('cfg must be a struct.');
end

moduleDir = fileparts(mfilename('fullpath'));
addpath(moduleDir);

repoRoot = fullfile(moduleDir, '..', '..', '..', '..');
cfg = fill_default_cfg(cfg, repoRoot);

if ~isfield(cfg, 'windowSec') || isempty(cfg.windowSec)
    error(['cfg.windowSec is required, for example cfg.windowSec = [4.32 8.64]. ', ...
           'The analysis window must be pre-specified or selected using training data only.']);
end

conditionMap = get_wm_condition_mapping(cfg.repoRoot);
cfg.conditionMap = conditionMap;
cfg.TR = detect_wm_tr(cfg.repoRoot, cfg);

[residCurves, subjectIDs, roiLabels, dataMetadata] = load_wm_roi_influence_data(cfg);
[nROI, nTime, nSubject, nCondition] = size(residCurves);
if nCondition ~= numel(conditionMap.conditionLabels)
    error('Expected %d WM conditions from repository mapping, found %d.', ...
        numel(conditionMap.conditionLabels), nCondition);
end

[windowTRs, windowSeconds] = window_seconds_to_trs(cfg.windowSec, cfg.TR, nTime);

pairedData = build_paired_interaction_data(residCurves, cfg, windowTRs, conditionMap);
if size(pairedData.D, 1) ~= nSubject || size(pairedData.D, 2) ~= nROI
    error('Paired data D has unexpected size [%s]; expected [%d x %d].', ...
        num2str(size(pairedData.D)), nSubject, nROI);
end

if ~exist(cfg.output.rootDir, 'dir')
    mkdir(cfg.output.rootDir);
end
splitFile = fullfile(cfg.output.rootDir, 'shared_subject_splits.mat');
outerSplits = make_repeated_subject_splits(nSubject, subjectIDs, cfg, splitFile);

runDir = fullfile(cfg.output.rootDir, sanitize_name(pairedData.contrastLabel), cfg.hrfModelName);
if ~exist(runDir, 'dir')
    mkdir(runDir);
end

checkpointFile = fullfile(runDir, 'roi_influence_checkpoint.mat');
outerResults = repmat(empty_outer_result(nROI, cfg.importance.nRandomSwaps), cfg.outer.nRepeats, 1);
completedRepeats = false(cfg.outer.nRepeats, 1);

if cfg.resume && isfile(checkpointFile)
    C = load(checkpointFile, 'outerResults', 'completedRepeats');
    if isfield(C, 'outerResults') && numel(C.outerResults) == cfg.outer.nRepeats
        outerResults = C.outerResults;
        completedRepeats = C.completedRepeats;
        fprintf('Resuming from checkpoint: %d/%d repeats complete.\n', ...
            sum(completedRepeats), cfg.outer.nRepeats);
    end
end

batchSize = cfg.outer.checkpointBatchSize;
for batchStart = 1:batchSize:cfg.outer.nRepeats
    batchEnd = min(cfg.outer.nRepeats, batchStart + batchSize - 1);
    batchIdx = find(~completedRepeats(batchStart:batchEnd)) + batchStart - 1;
    if isempty(batchIdx)
        continue
    end

    fprintf('Running outer repeats %d-%d (%d pending).\n', batchStart, batchEnd, numel(batchIdx));

    if cfg.outer.useParfor
        batchResults = cell(numel(batchIdx), 1);
        parfor ii = 1:numel(batchIdx)
            r = batchIdx(ii);
            batchResults{ii} = fit_outer_elasticnet_model( ...
                pairedData.D, subjectIDs, outerSplits(r), cfg);
        end
        for ii = 1:numel(batchIdx)
            r = batchIdx(ii);
            outerResults(r) = batchResults{ii};
            completedRepeats(r) = true;
        end
    else
        for ii = 1:numel(batchIdx)
            r = batchIdx(ii);
            outerResults(r) = fit_outer_elasticnet_model( ...
                pairedData.D, subjectIDs, outerSplits(r), cfg);
            completedRepeats(r) = true;
        end
    end

    save(checkpointFile, 'cfg', 'windowTRs', 'windowSeconds', ...
        'outerSplits', 'outerResults', 'completedRepeats', '-v7.3');
end

if ~all(completedRepeats)
    error('Not all outer repeats completed. Check checkpoint file: %s', checkpointFile);
end

signFlipResults = run_paired_signflip_maxT(pairedData.D, cfg);
summaryTable = summarize_roi_influence(outerResults, signFlipResults, roiLabels, cfg);

runMetadata = struct();
runMetadata.module = 'wm_roi_influence_elasticnet';
runMetadata.hrfModelName = cfg.hrfModelName;
runMetadata.dataFile = dataMetadata.dataFile;
runMetadata.dataVariable = dataMetadata.dataVariable;
runMetadata.dataSize = [nROI nTime nSubject nCondition];
runMetadata.subjectIDSource = dataMetadata.subjectIDSource;
runMetadata.roiLabelSource = dataMetadata.roiLabelSource;
runMetadata.conditionMapSource = conditionMap.sourceFiles;
runMetadata.contrastLabel = pairedData.contrastLabel;
runMetadata.analysisMode = cfg.analysisMode;
runMetadata.importanceComputation = dataMetadata.importanceComputation;
runMetadata.thresholdInterpretation = ['Thresholds define a reproducibility-based candidate set; ', ...
    'they are not causal or biologically definitive conclusions.'];

softwareVersionInfo = collect_software_version_info();
timestamp = datestr(now, 30);

save_roi_influence_outputs(cfg, runDir, windowTRs, windowSeconds, subjectIDs, ...
    roiLabels, outerSplits, outerResults, summaryTable, signFlipResults, ...
    runMetadata, softwareVersionInfo, timestamp);

if cfg.plots.makePlots
    plot_roi_influence_summary(summaryTable, outerResults, cfg, runDir);
end

out = struct();
out.cfg = cfg;
out.windowTRs = windowTRs;
out.windowSeconds = windowSeconds;
out.subjectIDs = subjectIDs;
out.roiLabels = roiLabels;
out.outerSplits = outerSplits;
out.outerResults = outerResults;
out.summaryTable = summaryTable;
out.signFlipResults = signFlipResults;
out.runMetadata = runMetadata;
out.softwareVersionInfo = softwareVersionInfo;
out.timestamp = timestamp;
out.runDir = runDir;

fprintf('ROI influence analysis complete: %s\n', runDir);
end

function cfg = fill_default_cfg(cfg, repoRoot)
cfg = set_default(cfg, 'repoRoot', repoRoot);
cfg = set_default(cfg, 'hrfModelName', 'cHRF');
cfg = set_default(cfg, 'contrastName', 'Body_LoadDiff_vs_Face_LoadDiff');
cfg = set_default(cfg, 'analysisMode', 'windowAverage');
cfg = set_default(cfg, 'seed', 20260621);
cfg = set_default(cfg, 'resume', true);

cfg = set_default(cfg, 'dataDir', fullfile(cfg.repoRoot, 'data', 'task_residual'));

if ~isfield(cfg, 'output') || ~isstruct(cfg.output)
    cfg.output = struct();
end
cfg.output = set_default(cfg.output, 'rootDir', ...
    fullfile(cfg.repoRoot, 'data', 'task_residual', 'roi_influence_elasticnet'));

if ~isfield(cfg, 'outer') || ~isstruct(cfg.outer)
    cfg.outer = struct();
end
cfg.outer = set_default(cfg.outer, 'nRepeats', 200);
cfg.outer = set_default(cfg.outer, 'trainFraction', 0.80);
cfg.outer = set_default(cfg.outer, 'checkpointBatchSize', 10);
cfg.outer = set_default(cfg.outer, 'useParfor', false);

if ~isfield(cfg, 'model') || ~isstruct(cfg.model)
    cfg.model = struct();
end
cfg.model = set_default(cfg.model, 'alpha', 0.50);
cfg.model = set_default(cfg.model, 'numLambda', 50);
cfg.model = set_default(cfg.model, 'innerKFolds', 5);
cfg.model = set_default(cfg.model, 'useOneSE', true);

if ~isfield(cfg, 'importance') || ~isstruct(cfg.importance)
    cfg.importance = struct();
end
cfg.importance = set_default(cfg.importance, 'nRandomSwaps', 50);
cfg.importance = set_default(cfg.importance, 'selectedOnly', true);

if ~isfield(cfg, 'permutation') || ~isstruct(cfg.permutation)
    cfg.permutation = struct();
end
cfg.permutation = set_default(cfg.permutation, 'nPerm', 5000);
cfg.permutation = set_default(cfg.permutation, 'alpha', 0.05);

if ~isfield(cfg, 'thresholds') || ~isstruct(cfg.thresholds)
    cfg.thresholds = struct();
end
cfg.thresholds = set_default(cfg.thresholds, 'minSelectionFrequency', 0.70);
cfg.thresholds = set_default(cfg.thresholds, 'minSignConsistency', 0.90);
cfg.thresholds = set_default(cfg.thresholds, 'minMeanReliance', 0);
cfg.thresholds = set_default(cfg.thresholds, 'requireFWER', true);

if ~isfield(cfg, 'plots') || ~isstruct(cfg.plots)
    cfg.plots = struct();
end
cfg.plots = set_default(cfg.plots, 'makePlots', true);
cfg.plots = set_default(cfg.plots, 'topN', 30);

if ~strcmpi(cfg.analysisMode, 'windowAverage')
    error('Only cfg.analysisMode = ''windowAverage'' is implemented. Time-resolved mode is reserved for later.');
end
end

function s = set_default(s, fieldName, value)
if ~isfield(s, fieldName) || isempty(s.(fieldName))
    s.(fieldName) = value;
end
end

function [windowTRs, windowSeconds] = window_seconds_to_trs(windowSec, TR, nTime)
if ~isnumeric(windowSec) || numel(windowSec) ~= 2 || any(~isfinite(windowSec))
    error('cfg.windowSec must be a finite numeric [startSec endSec] vector.');
end
windowSec = double(windowSec(:))';
if windowSec(2) < windowSec(1)
    error('cfg.windowSec endSec must be >= startSec.');
end
tSec = (0:(nTime-1)) * TR;
windowTRs = find(tSec >= windowSec(1) & tSec <= windowSec(2));
if isempty(windowTRs)
    error('cfg.windowSec [%g %g] contains no TRs for TR=%g and nTime=%d.', ...
        windowSec(1), windowSec(2), TR, nTime);
end
windowSeconds = tSec(windowTRs);
end

function e = empty_outer_result(nROI, nSwaps)
e = struct();
e.repeat = NaN;
e.seed = NaN;
e.betaRaw = nan(nROI, 1);
e.betaStandardized = nan(nROI, 1);
e.intercept = NaN;
e.selectedLambda = NaN;
e.trainSubjectIDs = [];
e.testSubjectIDs = [];
e.trainSubjectIdx = [];
e.testSubjectIdx = [];
e.TAFC = NaN;
e.testMargins = [];
e.selectionIndicator = false(nROI, 1);
e.haufePatternRaw = nan(nROI, 1);
e.haufePatternZ = nan(nROI, 1);
e.modelRelianceTAFC = zeros(nROI, nSwaps, 'single');
e.featureMedian = nan(1, nROI);
e.featureMean = nan(1, nROI);
e.featureScale = nan(1, nROI);
e.innerFoldIDByTrainSubject = [];
e.innerValidationLossMean = [];
e.innerValidationLossSE = [];
e.lambdaPath = [];
end

function s = sanitize_name(s)
s = char(s);
s = regexprep(s, '[^A-Za-z0-9_+-]+', '_');
s = regexprep(s, '_+', '_');
s = regexprep(s, '^_|_$', '');
end

function info = collect_software_version_info()
info = struct();
info.matlab = version;
info.computer = computer;
info.date = datestr(now, 30);
info.statistics = try_ver('stats');
info.parallel = try_ver('parallel');
end

function v = try_ver(toolboxName)
try
    v = ver(toolboxName);
catch
    v = [];
end
end
