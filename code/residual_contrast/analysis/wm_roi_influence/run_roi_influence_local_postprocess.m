function outputs = run_roi_influence_local_postprocess(resultFiles, varargin)
%RUN_ROI_INFLUENCE_LOCAL_POSTPROCESS Make light outputs from completed runs.
%
% Use this locally or on a login node after the cluster has produced the
% heavy roi_influence_*.mat file(s). This function does not refit models,
% rerun permutations, or touch checkpoints.
%
% Examples:
%   run_roi_influence_local_postprocess('/path/to/roi_influence_cHRF_*.mat');
%
%   files = struct();
%   files.cHRF = '/.../Body_LoadDiff_vs_Face_LoadDiff/cHRF/roi_influence_cHRF_Body_LoadDiff_vs_Face_LoadDiff.mat';
%   files.cHRFderiv = '/.../Body_LoadDiff_vs_Face_LoadDiff/cHRFderiv/roi_influence_cHRFderiv_Body_LoadDiff_vs_Face_LoadDiff.mat';
%   files.sHRF = '/.../Body_LoadDiff_vs_Face_LoadDiff/sHRF/roi_influence_sHRF_Body_LoadDiff_vs_Face_LoadDiff.mat';
%   run_roi_influence_local_postprocess(files);

p = inputParser;
addParameter(p, 'MakePlots', true, @islogical);
addParameter(p, 'RewriteCsv', true, @islogical);
addParameter(p, 'TopN', [], @(x) isempty(x) || (isscalar(x) && isnumeric(x) && x > 0));
addParameter(p, 'PlotSubdir', 'plots', @(x) ischar(x) || isstring(x));
addParameter(p, 'CompareModels', true, @islogical);
addParameter(p, 'ComparisonCsv', '', @(x) ischar(x) || isstring(x));
parse(p, varargin{:});
opt = p.Results;

files = normalize_result_files(resultFiles);
outputs = struct();
outputs.processedFiles = files(:);
outputs.plotDirs = cell(numel(files), 1);
outputs.csvFiles = cell(numel(files), 1);

for i = 1:numel(files)
    resultFile = files{i};
    S = load(resultFile, 'cfg', 'summaryTable', 'outerResults', 'runMetadata');
    required = {'cfg', 'summaryTable', 'outerResults'};
    for k = 1:numel(required)
        if ~isfield(S, required{k})
            error('Result file %s is missing %s.', resultFile, required{k});
        end
    end

    cfg = S.cfg;
    if ~isempty(opt.TopN)
        cfg.plots.topN = opt.TopN;
    elseif ~isfield(cfg, 'plots') || ~isfield(cfg.plots, 'topN') || isempty(cfg.plots.topN)
        cfg.plots.topN = 30;
    end

    runDir = fileparts(resultFile);
    plotDir = fullfile(runDir, char(opt.PlotSubdir));
    if opt.MakePlots
        plot_roi_influence_summary(S.summaryTable, S.outerResults, cfg, plotDir);
    end

    csvFile = fullfile(runDir, 'summaryTable.csv');
    if opt.RewriteCsv
        writetable(S.summaryTable, csvFile);
    end

    outputs.plotDirs{i} = plotDir;
    outputs.csvFiles{i} = csvFile;
end

outputs.comparisonCsv = '';
if opt.CompareModels && numel(files) == 3
    modelFiles = infer_model_file_struct(files);
    if isempty(char(opt.ComparisonCsv))
        contrastDir = fileparts(fileparts(modelFiles.cHRF));
        comparisonCsv = fullfile(contrastDir, 'hrf_model_comparison.csv');
    else
        comparisonCsv = char(opt.ComparisonCsv);
    end
    compare_roi_influence_hrf_models(modelFiles, comparisonCsv);
    outputs.comparisonCsv = comparisonCsv;
end

fprintf('Local ROI influence postprocess complete.\n');
end

function files = normalize_result_files(resultFiles)
if nargin < 1 || isempty(resultFiles)
    envFile = getenv('RESULT_MAT');
    if isempty(envFile)
        error('Provide resultFiles or set RESULT_MAT to a completed roi_influence_*.mat file.');
    end
    resultFiles = envFile;
end

if isstruct(resultFiles)
    names = {'cHRF', 'cHRFderiv', 'sHRF'};
    files = cell(1, 3);
    for i = 1:3
        if ~isfield(resultFiles, names{i})
            error('resultFiles struct must contain %s.', names{i});
        end
        files{i} = char(resultFiles.(names{i}));
    end
elseif ischar(resultFiles) || isstring(resultFiles)
    files = cellstr(resultFiles);
elseif iscell(resultFiles)
    files = resultFiles(:)';
else
    error('resultFiles must be a char/string path, cell array, or struct.');
end

for i = 1:numel(files)
    files{i} = char(files{i});
    if ~isfile(files{i})
        error('Result file not found: %s', files{i});
    end
end
end

function modelFiles = infer_model_file_struct(files)
modelFiles = struct();
for i = 1:numel(files)
    S = load(files{i}, 'cfg');
    if ~isfield(S, 'cfg') || ~isfield(S.cfg, 'hrfModelName')
        error('Cannot infer HRF model from %s.', files{i});
    end
    switch lower(char(S.cfg.hrfModelName))
        case 'chrf'
            modelFiles.cHRF = files{i};
        case 'chrfderiv'
            modelFiles.cHRFderiv = files{i};
        case 'shrf'
            modelFiles.sHRF = files{i};
        otherwise
            error('Unknown hrfModelName in %s: %s', files{i}, S.cfg.hrfModelName);
    end
end

required = {'cHRF', 'cHRFderiv', 'sHRF'};
for i = 1:numel(required)
    if ~isfield(modelFiles, required{i})
        error('Need cHRF, cHRFderiv, and sHRF result files for HRF comparison.');
    end
end
end
