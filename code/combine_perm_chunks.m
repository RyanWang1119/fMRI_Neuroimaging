function combine_perm_chunks(OUTDIR, TASKNAME, observed_mat_file, alpha)
if nargin < 4 || isempty(alpha), alpha = 0.05; end
TASKNAME = upper(TASKNAME);
pat = fullfile(OUTDIR, sprintf('perm_chunk_%s_chunk*.mat', lower(TASKNAME)));
d = dir(pat); assert(~isempty(d),'No permutation chunks found: %s', pat);

% load and stack maxima
maxTA  = []; maxTAFC = []; maxAbsH_global = []; maxAbsH_byRegion = [];
for i=1:numel(d)
    S = load(fullfile(d(i).folder, d(i).name));
    maxTA  = [maxTA;  S.maxTA_perm_chunk(:)];
    maxTAFC= [maxTAFC;S.maxTAFC_perm_chunk(:)];
    maxAbsH_global = [maxAbsH_global; S.maxAbsH_global_perm_chunk(:)];
    maxAbsH_byRegion = [maxAbsH_byRegion; S.maxAbsH_byRegion_perm_chunk]; %#ok<AGROW>
end
nPermTotal = numel(maxTA);
fprintf('[%s] Combined %d permutations from %d chunks.\n', TASKNAME, nPermTotal, numel(d));

% thresholds (FWER α) using max-stat distributions
thr.TA    = quantile(maxTA,   1-alpha);
thr.TAFC  = quantile(maxTAFC, 1-alpha);
thr.Hglob = quantile(maxAbsH_global, 1-alpha);
thr.HbyReg= quantile(maxAbsH_byRegion, 1-alpha, 1);   % 1 x Rg (per-region FWER across t×contrast)

% load observed to compute corrected significance & p-values
Obs = load(observed_mat_file, 'TA_mean','TAFC_mean','Haufe_mean','contrast_labels','cond_labels','T','Rg','S');
TA_obs   = Obs.TA_mean;     % [T x nCon]
TAFC_obs = Obs.TAFC_mean;   % [T x nCon]
H_mean   = Obs.Haufe_mean;  % [Rg x T x nCon]

% corrected significance masks
sig.TA   = TA_obs   >= thr.TA;
sig.TAFC = TAFC_obs >= thr.TAFC;

% global-map correction
H_abs = abs(H_mean);
sig.Haufe_global = H_abs >= thr.Hglob;

% by-region correction (across t×contrast only; controls FWER per region)
thr_vec = reshape(thr.HbyReg, [], 1); % [Rg x 1]
sig.Haufe_byRegion = bsxfun(@ge, H_abs, thr_vec);

% compute max-stat p-values for time/contrast cells
p.TA   = (1 + arrayfun(@(x) sum(maxTA   >= x), TA_obs))   ./ (1 + nPermTotal);
p.TAFC = (1 + arrayfun(@(x) sum(maxTAFC >= x), TAFC_obs)) ./ (1 + nPermTotal);

% save results
outf = fullfile(OUTDIR, sprintf('stageB_results_%s.mat', lower(TASKNAME)));
save(outf, 'thr','sig','p','nPermTotal', '-v7.3');
fprintf('[%s] Stage B combined saved: %s\n', TASKNAME, outf);
end
