addpath('/users/rwang');  
addpath(genpath('/users/rwang/CanLabCore'));           
addpath(genpath('/users/rwang/matlab_utilities'));    
addpath(genpath('/users/rwang/fMRI_Neuroimaging'));  
rehash path

% sanity on client: xval_SVM must be resolvable
assert(exist('xval_SVM','file')==2, 'xval_SVM.m not found on client path');

%% ---------- Robust parpool honoring SLURM + per-task JobStorageLocation ----------
w_env = str2double(getenv('SLURM_CPUS_PER_TASK'));
if isnan(w_env) || w_env < 1, w_env = feature('numcores'); end
maxNumCompThreads(1);

% Unique, writable job storage per array task
user   = getenv('USER'); if isempty(user), user = 'user'; end
jobid  = getenv('SLURM_JOB_ID');           if isempty(jobid),  jobid  = datestr(now,'yyyymmddHHMMSS'); end
taskid = getenv('SLURM_ARRAY_TASK_ID');    if isempty(taskid), taskid = '0'; end

baseJS = getenv('TMPDIR');                             % node-local scratch if available
if isempty(baseJS) || ~exist(baseJS,'dir')
    baseJS = fullfile('/users', user, 'tmp');          % fallback in $HOME
end
jobStorage = fullfile(baseJS, sprintf('matlab_jobstorage_%s_%s_%s', user, jobid, taskid));
if ~exist(jobStorage,'dir'); mkdir(jobStorage); end

pool = gcp('nocreate');
if isempty(pool)
    c = parcluster('local');
    try
        c.JobStorageLocation = jobStorage;             % isolate per task
    catch
        warning('Could not set JobStorageLocation; proceeding with default.');
    end
    allowed = c.NumWorkers;
    req = max(1, min([w_env, allowed]));
    try, c.NumWorkers = req; catch, req = min(req, c.NumWorkers); end

    fprintf('Starting parpool(local) with %d workers (SLURM=%d, cap=%d)\n', req, w_env, allowed);
    try
        pool = parpool(c, req);
    catch ME
        warning('parpool failed (%s). Retrying with fewer workers...', ME.message);
        req2 = max(1, floor(req/2));
        try, c.NumWorkers = req2; catch, req2 = min(req2, c.NumWorkers); end
        fprintf('Retrying parpool with %d workers.\n', req2);
        pool = parpool(c, req2);
    end
else
    fprintf('Using existing parpool with %d workers.\n', pool.NumWorkers);
end

% Propagate paths to all workers — so parfor can see xval_SVM, etc.
workerPathRoots = { ...
  '/users/rwang', ...
  '/users/rwang/CanLabCore', ...
  '/users/rwang/matlab_utilities', ...
  '/users/rwang/fMRI_Neuroimaging' ...
};
for i = 1:numel(workerPathRoots)
    root = workerPathRoots{i};
    if exist(root,'dir'), pctRunOnAll addpath(genpath(root)); end
end
pctRunOnAll rehash path
spmd, assert(exist('xval_SVM','file')==2, 'xval_SVM.m not found on worker'); end

%% ---------- Data file mappings (adjust paths if needed) ----------
files.EMOTION    = '/users/rwang/HRF_Resid_EMOTION_LR.mat';
files.GAMBLING   = '/users/rwang/HRF_Resid_GAMBLING_LR.mat';
files.LANGUAGE   = '/users/rwang/HRF_Resid_LANGUAGE_LR.mat';
files.MOTOR      = '/users/rwang/HRF_Resid_MOTOR_LR.mat';
files.RELATIONAL = '/users/rwang/HRF_Resid_RELATIONAL_LR.mat';
files.WM         = '/users/rwang/HRF_Resid_WMblock_LR.mat';

%% ---------- Outputs ----------
OUTDIR = getenv('OUTDIR');
if isempty(OUTDIR), OUTDIR = fullfile('/users', user, 'jhpce_outputs_stageA'); end
if ~exist(OUTDIR,'dir'), mkdir(OUTDIR); end

TASK = upper(string(getenv('TASK')));  % "" means run all
fprintf('TASK="%s"\nOUTDIR="%s"\nJobStorage="%s"\n', TASK, OUTDIR, jobStorage);

%% ---------- Dispatch ----------
runOne = @(T) dispatch_task(char(T), files, OUTDIR);

if TASK == ""
    order = ["EMOTION","GAMBLING","LANGUAGE","MOTOR","RELATIONAL","WM"];
    for t = order
        runOne(t);
    end
else
    runOne(TASK);
end

disp('All requested Stage A runs finished.');

%% ================== local function ===================
function dispatch_task(TASK, files, OUTDIR)
    fprintf('\n=== Running %s ===\n', TASK);
    switch TASK
        case 'EMOTION',    PredictContrast_EMOTION(files.EMOTION, OUTDIR);
        case 'GAMBLING',   PredictContrast_GAMBLING(files.GAMBLING, OUTDIR);
        case 'LANGUAGE',   PredictContrast_LANGUAGE(files.LANGUAGE, OUTDIR);
        case 'MOTOR',      PredictContrast_MOTOR(files.MOTOR, OUTDIR);
        case 'RELATIONAL', PredictContrast_RELATIONAL(files.RELATIONAL, OUTDIR);
        case 'WM',         PredictContrast_WM(files.WM, OUTDIR);
        otherwise, error('Unknown TASK=%s', TASK);
    end
end
