function out = run_roi_influence_elasticnet_slurm()
%RUN_ROI_INFLUENCE_ELASTICNET_SLURM Environment-driven wrapper for clusters.
%
% Environment variables:
%   HRF_MODEL       cHRF, cHRFderiv, or sHRF
%   CONTRAST_NAME   optional contrast label, e.g. Body_LoadDiff_vs_Face_LoadDiff
%   CONTRAST_INDEX  optional integer 1..6
%   WINDOW_SEC      required "start end", e.g. "4.32 8.64"
%   OUTDIR          optional output root
%   DATA_DIR        optional folder containing WMcHRF.mat, WMcHRFderiv.mat, WMsHRF.mat
%   NREPEATS        optional outer repeat count
%   NPERM           optional sign-flip permutation count
%   USE_PARFOR      optional 1/0

moduleDir = fileparts(mfilename('fullpath'));
addpath(moduleDir);

cfg = struct();
cfg.hrfModelName = getenv_default('HRF_MODEL', 'cHRF');

dataDir = getenv('DATA_DIR');
if ~isempty(dataDir)
    cfg.dataDir = dataDir;
end

contrastIndex = str2double(getenv('CONTRAST_INDEX'));
if isfinite(contrastIndex)
    cfg.contrastIndex = contrastIndex;
else
    cfg.contrastName = getenv_default('CONTRAST_NAME', 'Body_LoadDiff_vs_Face_LoadDiff');
end

windowText = getenv('WINDOW_SEC');
if isempty(windowText)
    error('WINDOW_SEC is required, for example WINDOW_SEC="4.32 8.64".');
end
cfg.windowSec = sscanf(windowText, '%f')';
if numel(cfg.windowSec) ~= 2
    error('WINDOW_SEC must contain two numeric values, e.g. "4.32 8.64".');
end

outdir = getenv('OUTDIR');
if ~isempty(outdir)
    cfg.output.rootDir = outdir;
end

nRepeats = str2double(getenv('NREPEATS'));
if isfinite(nRepeats) && nRepeats > 0
    cfg.outer.nRepeats = nRepeats;
end

nPerm = str2double(getenv('NPERM'));
if isfinite(nPerm) && nPerm > 0
    cfg.permutation.nPerm = nPerm;
end

useParfor = str2double(getenv('USE_PARFOR'));
if isfinite(useParfor)
    cfg.outer.useParfor = useParfor ~= 0;
end

if isfield(cfg, 'outer') && isfield(cfg.outer, 'useParfor') && cfg.outer.useParfor
    start_slurm_parpool();
end

out = run_roi_influence_elasticnet(cfg);
end

function val = getenv_default(name, defaultVal)
val = getenv(name);
if isempty(val)
    val = defaultVal;
end
end

function start_slurm_parpool()
pool = gcp('nocreate');
if ~isempty(pool)
    return
end

wEnv = str2double(getenv('SLURM_CPUS_PER_TASK'));
if ~isfinite(wEnv) || wEnv < 1
    wEnv = feature('numcores');
end

c = parcluster('local');
req = max(1, min(wEnv, c.NumWorkers));

user = getenv('USER');
if isempty(user)
    user = getenv('USERNAME');
end
if isempty(user)
    user = 'user';
end

jobid = getenv('SLURM_JOB_ID');
if isempty(jobid)
    jobid = datestr(now, 'yyyymmddHHMMSS');
end
taskid = getenv('SLURM_ARRAY_TASK_ID');
if isempty(taskid)
    taskid = '0';
end
baseJS = getenv('TMPDIR');
if isempty(baseJS) || ~exist(baseJS, 'dir')
    baseJS = tempdir;
end
jobStorage = fullfile(baseJS, sprintf('matlab_roi_influence_%s_%s_%s', user, jobid, taskid));
if ~exist(jobStorage, 'dir')
    mkdir(jobStorage);
end

try
    c.JobStorageLocation = jobStorage;
catch
    warning('Could not set JobStorageLocation; using MATLAB default.');
end

fprintf('Starting parpool(local) with %d workers.\n', req);
parpool(c, req);
end
