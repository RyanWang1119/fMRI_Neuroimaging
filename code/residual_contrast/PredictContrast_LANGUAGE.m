function PredictContrast_LANGUAGE(DATA_FILE, OUTDIR, OUTTAG)

if nargin < 2 || isempty(OUTDIR), OUTDIR = fullfile(pwd,'outputs'); end
if nargin < 3, OUTTAG = ''; end
if ~exist(OUTDIR,'dir'), mkdir(OUTDIR); end

% ---- Load & robust remap to [S x R x T x C] ----
tmp = load(DATA_FILE, 'Results');
assert(isfield(tmp, 'Results'), 'Results missing in %s', DATA_FILE);
A = tmp.Results; clear tmp

sz = size(A); sz(end+1:4) = 1;         
dims = 1:4;

iR = find(sz == 489, 1, 'first');
iT = find(sz == 41 , 1, 'first');
iC = find(sz == 6  , 1, 'first');
assert(~isempty(iR) && ~isempty(iT) && ~isempty(iC), ...
    'LANG dims not recognizable: size=%s (need one 489, one 41, one 6).', mat2str(sz));

iS = setdiff(dims, [iR iT iC]);         
assert(numel(iS) == 1, 'Could not infer a unique subjects axis from size=%s', mat2str(sz));

% Permute to [S R T C]
R = permute(A, [iS iR iT iC]);           % [S x 489 x 41 x 6]
[S, Rg, T, C] = size(R);
assert(Rg == 489 && T == 41 && C == 6, ...
    'After remap expected [S x 489 x 41 x 6], got [%d x %d x %d x %d].', S, Rg, T, C);

% ---- Construct math–story diffs within P/Q/R ----
% HCP order: 1..3 = math [present, question, response]; 4..6 = story [present, question, response]
Diff_P = R(:,:,:,1) - R(:,:,:,4);   % present:  math - story
Diff_Q = R(:,:,:,2) - R(:,:,:,5);   % question: math - story
Diff_R = R(:,:,:,3) - R(:,:,:,6);   % response: math - story

Results = cat(4, Diff_P, Diff_Q, Diff_R);  % [S x 489 x 41 x 3]
cond_labels = {'Present_MathMinusStory','Question_MathMinusStory','Response_MathMinusStory'};
COMP2 = nchoosek(1:3, 2);                 % pairwise among the 3 diffs

stamp = datestr(now, 'HHMMSSFFF');
if ~isempty(OUTTAG)
    temp_file = fullfile(OUTDIR, sprintf('temp_LANG_%s_%s.mat', stamp, OUTTAG));
    taskname  = ['LANGUAGE_Interaction_' OUTTAG];
else
    temp_file = fullfile(OUTDIR, sprintf('temp_LANG_%s.mat', stamp));
    taskname  = 'LANGUAGE_Interaction';
end

AlreadyCanonical = true; %#ok<NASGU>  
save(temp_file, 'Results', 'AlreadyCanonical', '-v7.3');

try
    run_observed_stage_flexible(temp_file, OUTDIR, COMP2, cond_labels, taskname);
catch ME
    if exist(temp_file, 'file'), delete(temp_file); end
    rethrow(ME);
end

if exist(temp_file, 'file'), delete(temp_file); end
end
