clear; clc;

%% ------------------ Config & paths --------------------------------------
STAGE     = getenv_default("STAGE","observed");  % observed | perm | merge
DATA_FILE = getenv_default("DATA_FILE", fullfile(pwd,'HRF_Resid_WMblock_LR.mat'));
OUTDIR    = getenv_default("OUTDIR", fullfile(getenv_default("HOME",pwd), 'jhpce_outputs_wmblock'));
if ~exist(OUTDIR,'dir'), mkdir(OUTDIR); end

fprintf('Stage      : %s\n', STAGE);
fprintf('Data file  : %s\n', DATA_FILE);
fprintf('Output dir : %s\n', OUTDIR);

% Worker setup from SLURM; avoid thread oversubscription
w = str2double(getenv('SLURM_CPUS_PER_TASK')); if isnan(w) || w < 1, w = 24; end
maxNumCompThreads(1);
if isempty(gcp('nocreate')), parpool('Processes', w); end
fprintf('Using %d MATLAB workers\n', w);

%% ------------------ Constants / design ----------------------------------
regions     = 489;
time_points = 41;
subjects    = 393;

CONDMAT = [5 1 6 2;   % body vs faces
           5 1 7 3;   % body vs places
           5 1 8 4;   % body vs tools
           6 2 7 3;   % faces vs places
           6 2 8 4;   % faces vs tools
           7 3 8 4];  % places vs tools
nCon = size(CONDMAT,1);

% Core controls
nBoot  = 20;      % bootstrap replicates
trainN = 314;     % training subjects
testN  = subjects - trainN;

%% ------------------ Stage switch ----------------------------------------
switch lower(STAGE)
    case 'observed'
        run_observed_stage(DATA_FILE, OUTDIR, regions, time_points, subjects, nCon, trainN, testN, nBoot, CONDMAT); 
    case 'perm'
        run_permutation_chunk(DATA_FILE, OUTDIR, regions, time_points, subjects, nCon); 
    case 'merge'
        run_merge_stage(OUTDIR);
    otherwise
        error('Unknown STAGE=%s (use observed|perm|merge)', STAGE);
end

disp('Done.');

%% ========================================================================
%                              STAGE A
% ========================================================================
function run_observed_stage(DATA_FILE, OUTDIR, regions, time_points, subjects, nCon, trainN, testN, nBoot, CONDMAT) %#ok<INUSD>
    % Precompute splits & bootstrap indices on client (saved for perm/merge)
    rng(42,'twister');
    trainSubs = cell(nBoot,1);
    testSubs  = cell(nBoot,1);
    for b = 1:nBoot
        trainSubs{b} = 1:trainN;
        testSubs{b}  = trainN+1:subjects;
    end
    bs_idx_list = cell(nBoot,1);
    for b = 1:nBoot
        bs_idx_list{b} = randsample(trainN, trainN, true);
    end

    Rconst = parallel.pool.Constant( ...
        @() load(getenv_default("DATA_FILE", fullfile(pwd,'HRF_Resid_WMblock_LR.mat')), 'Results'), ...
        @(S) permute(S.Results, [3 1 2 4]) ...
    );
    CONDMATconst   = parallel.pool.Constant(CONDMAT);
    trainSubsConst = parallel.pool.Constant(trainSubs);
    testSubsConst  = parallel.pool.Constant(testSubs);
    bsIdxConst     = parallel.pool.Constant(bs_idx_list);

    PERM_TOTAL = str2double(getenv_default("PERM_TOTAL","10000"));

    % ================== Preallocate observed outputs ======================
    % TA/TAFC in both layouts (time x contrast x boot) and (time x boot x contrast)
    TA_obs      = zeros(time_points, nCon, nBoot);
    TAFC_obs    = zeros(time_points, nCon, nBoot);
    TA_boot     = zeros(time_points, nBoot, nCon);   
    TAFC_boot   = zeros(time_points, nBoot, nCon);   

    % Weights: region-major (legacy) and Beta_obs (time x boot x contrast x region)
    W_obs       = zeros(regions, time_points, nCon, nBoot, 'single');     % [reg x time x con x boot]
    Haufe_obs   = zeros(regions, time_points, nCon, nBoot, 'single');
    Beta_obs    = zeros(time_points, nBoot, nCon, regions, 'single');     % [time x boot x con x reg]

    parfor b = 1:nBoot
        rng(1000+b,'twister');

        % Pull worker-local constants
        R        = Rconst.Value;
        CONDMAT_ = CONDMATconst.Value;
        tr       = trainSubsConst.Value{b};
        te       = testSubsConst.Value{b};
        bs_idx   = bsIdxConst.Value{b};

        TA_b    = zeros(time_points, nCon);      % [time x con]
        TAFC_b  = zeros(time_points, nCon);      % [time x con]
        W_b     = zeros(regions, time_points, nCon, 'single'); % [reg x time x con]
        H_b     = zeros(regions, time_points, nCon, 'single'); % [reg x time x con]

        for c = 1:nCon
            A = CONDMAT_(c,1); B = CONDMAT_(c,2);
            C = CONDMAT_(c,3); D = CONDMAT_(c,4);

            X1_all = squeeze(R(:,:,:,A) - R(:,:,:,B)); % [sub x reg x time]
            X2_all = squeeze(R(:,:,:,C) - R(:,:,:,D));

            for t = 1:time_points
                X1t = X1_all(:,:,t); X2t = X2_all(:,:,t);

                % Train & test (paired)
                X1_tr = X1t(tr,:); X2_tr = X2t(tr,:);
                Xtr   = [X1_tr(bs_idx,:); X2_tr(bs_idx,:)];
                Ytr   = [ ones(numel(bs_idx),1); -ones(numel(bs_idx),1) ];
                tr_ids = [ tr(bs_idx)'; tr(bs_idx)' ]; %#ok<NASGU>

                Xte = [X1t(te,:); X2t(te,:)];
                Yte = [ ones(numel(te),1); -ones(numel(te),1) ];

                % Train-only preprocess (impute + drop zero-var + z-score)
                [XtrZ, mu_full, sd_full, keepIdx] = standardize_impute_train(Xtr);
                XteZ = standardize_impute_apply(Xte, mu_full, sd_full, keepIdx);

                % Train (xval_SVM) and weights
                Svm = xval_SVM(XtrZ, Ytr, [ tr(bs_idx)'; tr(bs_idx)' ], ...
                               'nooptimize','norepeats','nobootstrap','noverbose','noplot');

                w_keep = Svm.ClassificationModel.Beta;

                % Expand to full region vector in *region-major* layout (legacy)
                w_full = zeros(regions,1,'single'); w_full(keepIdx) = single(w_keep);
                W_b(:,t,c) = w_full;

                % Haufe (Z-space) -> back to original via sd
                CovZ   = cov(XtrZ);
                a_keep = CovZ * w_keep;
                a_full = zeros(regions,1,'single');
                a_full(keepIdx) = single( a_keep ./ sd_full(keepIdx)' );
                H_b(:,t,c) = a_full;

                % Predict / metrics
                [yhat, score] = predict(Svm.ClassificationModel, XteZ);
                TA_b(t,c) = mean(yhat == Yte);

                poscol = find(Svm.ClassificationModel.ClassNames==1,1);
                score_pos = score(:,poscol);
                s1 = score_pos(1:numel(te)); s2 = score_pos(numel(te)+1:end);
                TAFC_b(t,c) = mean(s1 > s2);
            end
        end

        % -------- parfor-safe slice assignments (match shapes exactly) -----
        % Time x contrast views (legacy)
        TA_obs(:,:,b)   = TA_b;                      % [time x con] -> [time x con x 1]
        TAFC_obs(:,:,b) = TAFC_b;
   
        TA_boot(:, b, :)   = reshape(TA_b,   [time_points, 1, nCon]);
        TAFC_boot(:, b, :) = reshape(TAFC_b, [time_points, 1, nCon]);

        % Weights/Haufe in region-major 
        W_obs(:,:,:,b)     = W_b;                    % [reg x time x con] -> ... x boot
        Haufe_obs(:,:,:,b) = H_b;

        % Also store weights as Beta in (time x boot x con x region)
        %   Convert W_b [reg x time x con] -> [time x con x reg], then slice at boot
        Beta_obs(:, b, :, :) = permute(W_b, [2 3 1]);
    end

    % Touch outputs (silence analyzer)
    tmpW = size(W_obs); %#ok<NASGU>
    tmpB = size(Beta_obs); %#ok<NASGU>

    % Aggregates (over boot)
    TA_mean   = mean(TA_obs,   3);
    TA_se     = std( TA_obs,0, 3) / sqrt(nBoot);
    TAFC_mean = mean(TAFC_obs, 3);
    TAFC_se   = std( TAFC_obs,0, 3) / sqrt(nBoot);
    Haufe_mean = mean(Haufe_obs, 4, 'omitnan');
    Haufe_std  = std( Haufe_obs, 0, 4, 'omitnan');

    % Save observed + metadata (no permSwap saved)
    obs_file = fullfile(OUTDIR, 'observed_wmblock.mat');
    save(obs_file, ...
        'TA_obs','TAFC_obs','TA_mean','TA_se','TAFC_mean','TAFC_se', ...
        'TA_boot','TAFC_boot','W_obs','Haufe_obs','Haufe_mean','Haufe_std', ...
        'Beta_obs', ... % NEW 4D (time x boot x con x region)
        'trainSubs','testSubs','bs_idx_list', ...
        'regions','time_points','subjects','nCon','trainN','testN','nBoot','CONDMAT', ...
        'PERM_TOTAL', ...
        '-v7.3');
    fprintf('Observed stage saved: %s\n', obs_file);
end

%% ========================================================================
%                              STAGE B
% ========================================================================
function run_permutation_chunk(DATA_FILE, OUTDIR, regions, time_points, subjects, nCon) %#ok<INUSD>
    % Load observed metadata (splits, bootstraps, design)
    obs_file = fullfile(OUTDIR, 'observed_wmblock.mat');
    Sobs = load(obs_file, ...
        'trainSubs','testSubs','bs_idx_list', ...
        'regions','time_points','subjects','nCon','nBoot','CONDMAT','PERM_TOTAL');
    trainSubs = Sobs.trainSubs; testSubs = Sobs.testSubs; bs_idx_list = Sobs.bs_idx_list;
    regions   = Sobs.regions;   time_points = Sobs.time_points;
    nCon      = Sobs.nCon;      nBoot = Sobs.nBoot; CONDMAT = Sobs.CONDMAT;
    PERM_TOTAL = Sobs.PERM_TOTAL;
    clear Sobs;

    % Workers load data directly (avoid broadcast)
    Rconst = parallel.pool.Constant( ...
        @() load(getenv_default("DATA_FILE", fullfile(pwd,'HRF_Resid_WMblock_LR.mat')), 'Results'), ...
        @(S) permute(S.Results, [3 1 2 4]) ...
    );
    CONDMATconst   = parallel.pool.Constant(CONDMAT);
    trainSubsConst = parallel.pool.Constant(trainSubs);
    testSubsConst  = parallel.pool.Constant(testSubs);
    bsIdxConst     = parallel.pool.Constant(bs_idx_list);

    % Perm chunk
    PERM_CHUNK = str2double(getenv_default("PERM_CHUNK","200"));
    taskIdx    = str2double(getenv_default("SLURM_ARRAY_TASK_ID","1"));
    permStart  = (taskIdx-1)*PERM_CHUNK + 1;
    permEnd    = min(permStart + PERM_CHUNK - 1, PERM_TOTAL);
    fprintf('Permutation chunk: %d..%d (of %d)\n', permStart, permEnd, PERM_TOTAL);
    nPermChunk = permEnd - permStart + 1;

    % Allocate outputs for this chunk
    maxTA_perm_chunk            = zeros(nPermChunk, nCon, nBoot);
    maxTAFC_perm_chunk          = zeros(nPermChunk, nCon, nBoot);
    maxAbsH_global_perm_chunk   = zeros(nPermChunk, nCon, nBoot);
    maxAbsH_byRegion_perm_chunk = zeros(regions, nPermChunk, nCon, nBoot, 'single');

    parfor pAbs = 1:nPermChunk
        p = permStart + pAbs - 1;
        rng(2000 + p, 'twister');

        % Pull worker-local constants
        R        = Rconst.Value;
        CONDMAT_ = CONDMATconst.Value;

        maxTA_p    = zeros(nCon, nBoot);
        maxTAFC_p  = zeros(nCon, nBoot);
        maxHglob_p = zeros(nCon, nBoot);
        maxHreg_p  = zeros(regions, nCon, nBoot, 'single');

        for b = 1:nBoot
            tr = trainSubsConst.Value{b};
            te = testSubsConst.Value{b};
            bs_idx = bsIdxConst.Value{b};

            for c = 1:nCon
                A = CONDMAT_(c,1); B = CONDMAT_(c,2);
                C = CONDMAT_(c,3); D = CONDMAT_(c,4);

                % Deterministic per-(b,c,p) swap (no permSwap broadcast)
                permSeed = 1e7 * b + 1e4 * c + p;
                rng(permSeed, 'twister');
                swap = rand(subjects,1) > 0.5;

                X1_all = squeeze(R(:,:,:,A) - R(:,:,:,B));
                X2_all = squeeze(R(:,:,:,C) - R(:,:,:,D));

                TA_t    = zeros(time_points,1);
                TAFC_t  = zeros(time_points,1);
                maxHreg = -inf(regions,1,'single');
                maxHglob= -inf;

                for t = 1:time_points
                    X1t = X1_all(:,:,t); X2t = X2_all(:,:,t);

                    % Within-subject swap
                    [X1_sw, X2_sw] = apply_swaps_pairwise(X1t, X2t, swap);

                    % Train & test (paired)
                    X1_tr = X1_sw(tr,:); X2_tr = X2_sw(tr,:);
                    Xtr   = [X1_tr(bs_idx,:); X2_tr(bs_idx,:)];
                    Ytr   = [ ones(numel(bs_idx),1); -ones(numel(bs_idx),1) ];
                    tr_ids = [ tr(bs_idx)'; tr(bs_idx)' ]; %#ok<NASGU>

                    Xte = [X1_sw(te,:); X2_sw(te,:)];
                    Yte = [ ones(numel(te),1); -ones(numel(te),1) ];

                    % Train-only preprocess
                    [XtrZ, mu_full, sd_full, keepIdx] = standardize_impute_train(Xtr);
                    XteZ = standardize_impute_apply(Xte, mu_full, sd_full, keepIdx);

                    % Train/eval
                    Svm = xval_SVM(XtrZ, Ytr, [ tr(bs_idx)'; tr(bs_idx)' ], ...
                                   'nooptimize','norepeats','nobootstrap','noverbose','noplot');
                    w_keep = Svm.ClassificationModel.Beta;

                    CovZ   = cov(XtrZ);
                    a_keep = CovZ * w_keep;
                    a_full = zeros(regions,1,'single');
                    a_full(keepIdx) = single( a_keep ./ sd_full(keepIdx)' );

                    [yhat, score] = predict(Svm.ClassificationModel, XteZ);
                    TA_t(t) = mean(yhat == Yte);

                    poscol = find(Svm.ClassificationModel.ClassNames==1,1);
                    score_pos = score(:,poscol);
                    s1 = score_pos(1:numel(te));
                    s2 = score_pos(numel(te)+1:end);
                    TAFC_t(t) = mean(s1 > s2);

                    % Haufe max trackers
                    absA = abs(a_full);
                    maxHreg = max(maxHreg, absA);
                    maxHglob = max(maxHglob, max(absA));
                end

                maxTA_p(c,b)    = max(TA_t);
                maxTAFC_p(c,b)  = max(TAFC_t);
                maxHglob_p(c,b) = maxHglob;
                maxHreg_p(:,c,b)= maxHreg;
            end
        end

        maxTA_perm_chunk(pAbs,:,:)              = maxTA_p;
        maxTAFC_perm_chunk(pAbs,:,:)            = maxTAFC_p;
        maxAbsH_global_perm_chunk(pAbs,:,:)     = maxHglob_p;
        maxAbsH_byRegion_perm_chunk(:,pAbs,:,:) = maxHreg_p;
    end

    % Touch outputs (silence analyzer)
    tmpP = size(maxTA_perm_chunk); %#ok<NASGU>

    % Save chunk
    chunk_file = fullfile(OUTDIR, sprintf('perm_chunk_%06d_%06d.mat', permStart, permEnd));
    save(chunk_file, ...
        'maxTA_perm_chunk','maxTAFC_perm_chunk', ...
        'maxAbsH_global_perm_chunk','maxAbsH_byRegion_perm_chunk', ...
        'permStart','permEnd','-v7.3');
    fprintf('Saved chunk: %s\n', chunk_file);
end

%% ========================================================================
%                              STAGE C
% ========================================================================
function run_merge_stage(OUTDIR)

    % ---- Load observed stats & metadata ----
    obs_file = fullfile(OUTDIR, 'observed_wmblock.mat');
    S = load(obs_file, ...
        'TA_obs','TAFC_obs','Haufe_obs', ...
        'regions','time_points','nCon','nBoot','PERM_TOTAL');
    TA_obs    = S.TA_obs;       % [time x con x boot]
    TAFC_obs  = S.TAFC_obs;     % [time x con x boot]
    Haufe_obs = S.Haufe_obs;    % [reg x time x con x boot], single
    regions     = S.regions;
    time_points = S.time_points;
    nCon        = S.nCon;
    nBoot       = S.nBoot;
    PERM_TOTAL  = S.PERM_TOTAL;

    % ---- Find permutation chunk files ----
    d = dir(fullfile(OUTDIR, 'perm_chunk_*.mat'));
    if isempty(d)
        error('No permutation chunks found in %s', OUTDIR);
    end
    [~,ord] = sort({d.name});
    d = d(ord);

    % ---- Allocate streaming counters (uint32) ----
    % TA/TAFC FWER across time: p per (t,c,b)
    cnt_TA   = zeros(time_points, nCon, nBoot, 'uint32');
    cnt_TAFC = zeros(time_points, nCon, nBoot, 'uint32');

    % Haufe: global FWER across time×regions (null is global max over time×regions)
    cnt_Hglob = zeros(regions, time_points, nCon, nBoot, 'uint32');

    % Haufe: per-region timewise FWER (null is per-region max over time)
    cnt_Hreg = zeros(regions, time_points, nCon, nBoot, 'uint32');

    % Track how many permutations we actually processed (robust to partial sets)
    nPerm_accum = 0;

    % ---- Stream over chunks ----
    for i = 1:numel(d)
        C = load(fullfile(OUTDIR,d(i).name), ...
            'maxTA_perm_chunk','maxTAFC_perm_chunk', ...
            'maxAbsH_global_perm_chunk','maxAbsH_byRegion_perm_chunk', ...
            'permStart','permEnd');

        maxTA_chunk     = C.maxTA_perm_chunk;           % [nPermChunk x nCon x nBoot]
        maxTAFC_chunk   = C.maxTAFC_perm_chunk;         % [nPermChunk x nCon x nBoot]
        maxHglob_chunk  = C.maxAbsH_global_perm_chunk;  % [nPermChunk x nCon x nBoot]
        maxHreg_chunk   = C.maxAbsH_byRegion_perm_chunk;% [regions x nPermChunk x nCon x nBoot], single
        nPermChunk      = size(maxTA_chunk, 1);

        % ---------- TA / TAFC: count exceedances per (t,c,b) ----------
        for b = 1:nBoot
          for c = 1:nCon
              nullTA   = squeeze(maxTA_chunk(:,  c, b));   % [nPermChunk x 1]
              nullTAFC = squeeze(maxTAFC_chunk(:,c, b));   % [nPermChunk x 1]
              % Vectorized across time points by simple loop (41*200 is tiny)
              for t = 1:time_points
                  cnt_TA(t,c,b)   = cnt_TA(t,c,b)   + uint32(sum(nullTA   >= TA_obs(t,c,b)));
                  cnt_TAFC(t,c,b) = cnt_TAFC(t,c,b) + uint32(sum(nullTAFC >= TAFC_obs(t,c,b)));
              end
          end
        end

        % ---------- Haufe global FWER: null is scalar max over time×regions ----------
        % For each (c,b), compare each scalar null against abs(H_obs)(:,t)
        for b = 1:nBoot
          for c = 1:nCon
              nullGlob = squeeze(maxHglob_chunk(:, c, b));   % [nPermChunk x 1], single/double
              % abs observed maps for this (c,b): [regions x time]
              absObs_rt = abs( single( Haufe_obs(:,:,c,b) ) );  % ensure single to save memory
              % Vectorized comparison: (nPermChunk x 1) >= (regions*time)
              % Build [nPermChunk x (regions*time)] by implicit expansion, then sum along rows.
              cmp_counts = sum( nullGlob >= reshape(absObs_rt, 1, []), 1 ); % 1 x (R*T)
              cnt_Hglob(:,:,c,b) = cnt_Hglob(:,:,c,b) + ...
                                   reshape( uint32(cmp_counts), regions, time_points );
          end
        end

        % ---------- Haufe per-region timewise FWER ----------
        % For each (c,b), compare region-wise null (regions x nPermChunk) with absObs (regions x time)
        for b = 1:nBoot
          for c = 1:nCon
              nullReg = single( squeeze(maxHreg_chunk(:,:,c,b)) ); % [regions x nPermChunk], single
              absObs_rt = abs( single( Haufe_obs(:,:,c,b) ) );     % [regions x time]
              cmp = bsxfun(@ge, nullReg, permute(absObs_rt, [1 3 2])); % [regions x nPermChunk x time]
              counts = squeeze( sum(cmp, 2) );                        % [regions x time]
              cnt_Hreg(:,:,c,b) = cnt_Hreg(:,:,c,b) + uint32(counts);
          end
        end

        nPerm_accum = nPerm_accum + nPermChunk;
        fprintf('Merged chunk %d/%d: perms so far = %d\n', i, numel(d), nPerm_accum);
    end

    % ---- Convert counts to p-values (unbiased +1 correction) ----
    denom = double(nPerm_accum) + 1.0;
    p_TA_FWER   = (double(cnt_TA)   + 1.0) / denom;   % [time x con x boot]
    p_TAFC_FWER = (double(cnt_TAFC) + 1.0) / denom;

    p_Haufe_FWER_global   = (double(cnt_Hglob) + 1.0) / denom; % [reg x time x con x boot]
    p_Haufe_FWER_byRegion = (double(cnt_Hreg)  + 1.0) / denom;

    % ---- Optional BH–FDR across regions at each (t,c,b) ----
    q = 0.05;
    Haufe_FDR_mask = false(regions, time_points, nCon, nBoot);
    for b = 1:nBoot
      for c = 1:nCon
        for t = 1:time_points
            pvec = p_Haufe_FWER_byRegion(:, t, c, b);
            [h,~,~] = fdr_bh(pvec, q, 'pdep', 'no'); % BH independence/PRDS
            Haufe_FDR_mask(:, t, c, b) = h;
        end
      end
    end

    % ---- Save final outputs ----
    final_file = fullfile(OUTDIR,'WMblock_residual_decoding_full.mat');
    save(final_file, ...
        'TA_obs','TAFC_obs', ...
        'p_TA_FWER','p_TAFC_FWER', ...
        'Haufe_obs','p_Haufe_FWER_global','p_Haufe_FWER_byRegion','Haufe_FDR_mask', ...
        'nPerm_accum', '-v7.3');
    fprintf('Merged results saved: %s (nPerm=%d)\n', final_file, nPerm_accum);
end


%% ========================================================================
%                             HELPERS
% ========================================================================
function v = getenv_default(name, defaultVal)
    v = getenv(name);
    if isempty(v), v = defaultVal; end
end

function [XtrZ, mu_full, sd_full, keepIdx] = standardize_impute_train(Xtr)
% Train-only: mean-impute NaNs, drop all-NaN/zero-var features, z-score.
    mu_full = mean(Xtr, 1, 'omitnan');                 % replaces nanmean
    Xtr_filled = Xtr;
    nanmask = isnan(Xtr_filled);
    if any(nanmask(:))
        for j = 1:size(Xtr_filled,2)
            if any(nanmask(:,j))
                Xtr_filled(nanmask(:,j), j) = mu_full(j);
            end
        end
    end
    sd_full = std(Xtr_filled, 0, 1);
    keepIdx = ~isnan(mu_full) & (sd_full > 0);
    if ~any(keepIdx)
        error('All features were NaN or zero-variance in training after imputation.');
    end
    mu_full(~keepIdx) = 0;
    sd_full(~keepIdx) = 1;
    XtrZ = (Xtr_filled(:, keepIdx) - mu_full(keepIdx)) ./ sd_full(keepIdx);
end

function XteZ = standardize_impute_apply(Xte, mu_full, sd_full, keepIdx)
% Apply the SAME imputation mask and standardization to test.
    Xte_filled = Xte;
    nanmask = isnan(Xte_filled);
    if any(nanmask(:))
        for j = 1:size(Xte_filled,2)
            if any(nanmask(:,j))
                Xte_filled(nanmask(:,j), j) = mu_full(j);
            end
        end
    end
    XteZ = (Xte_filled(:, keepIdx) - mu_full(keepIdx)) ./ sd_full(keepIdx);
end

function [X1_sw, X2_sw] = apply_swaps_pairwise(X1, X2, swapMask)
% Swap paired rows (subjects where swapMask==true).
    X1_sw = X1; X2_sw = X2;
    idx = find(swapMask);
    if ~isempty(idx)
        tmp = X1_sw(idx, :);
        X1_sw(idx, :) = X2_sw(idx, :);
        X2_sw(idx, :) = tmp;
    end
end

function [h, crit_p, adj_p] = fdr_bh(pvals, q, method, report)
% Benjamini-Hochberg FDR (indep/posdep)
    if nargin < 2 || isempty(q), q = 0.05; end
    if nargin < 3 || isempty(method), method = 'pdep'; end
    if nargin < 4, report = 'no'; end
    p = pvals(:);
    [p_sorted, sort_ids] = sort(p);
    m = numel(p);
    if strcmpi(method,'pdep')
        denom = 1:m;
    else
        denom = sum(1./(1:m)) * ones(1,m);
    end
    thresh = (1:m)' ./ denom' * q;
    below = p_sorted <= thresh;
    if any(below)
        kmax = find(below,1,'last');
        crit_p = p_sorted(kmax);
        h = p <= crit_p;
    else
        crit_p = 0;
        h = false(size(p));
    end
    % adjusted p (monotone)
    adj_p_sorted = zeros(size(p_sorted));
    for k = m:-1:1
        adj_p_sorted(k) = min( (m/denom(min(k,end))) * p_sorted(k), ...
                               (k<m)*adj_p_sorted(k+1) + (k==m)*1 );
    end
    adj_p = zeros(size(p));
    adj_p(sort_ids) = adj_p_sorted;
    if strcmpi(report,'yes')
        fprintf('FDR (BH-%s): %d/%d < q=%.3f\n', method, sum(h), m, q);
    end
end
