function interpret_observed_alltasks(taskInput, varargin)
% Loads observed_<task>.mat (Stage A) and produces:
%  1. Timecourse plots with Error Ribbons (TA & PMZ)
%  2. Summary CSVs of peak timing and magnitude
%  3. Spatial Stats CSV (Weights vs Haufe Z-scores)
%  4. Haufe Heatmaps for Top 50 regions

%% ---- USER SETTINGS ----
INDIR  = fullfile(pwd, 'data', 'observed_residual_results');
OUTDIR = fullfile(INDIR, 'observed_interpretation_2');
TR     = 0.72;

% Defaults (used only if taskInput not provided)
defaultTaskFiles = struct( ...
  'EMOTION',    fullfile(INDIR, 'observed_emotion.mat'), ...
  'GAMBLING',   fullfile(INDIR, 'observed_gambling.mat'), ...
  'LANGUAGE',   fullfile(INDIR, 'observed_language_interaction.mat'), ...
  'MOTOR',      fullfile(INDIR, 'observed_motor.mat'), ...
  'RELATIONAL', fullfile(INDIR, 'observed_relational.mat'), ...
  'WM',         fullfile(INDIR, 'observed_wm_interaction.mat') ...
);

% --- Parse optional name-value overrides ---
p = inputParser;
addParameter(p, 'INDIR',  INDIR);
addParameter(p, 'OUTDIR', OUTDIR);
addParameter(p, 'TR',     TR);
parse(p, varargin{:});
INDIR  = p.Results.INDIR;
OUTDIR = p.Results.OUTDIR;
TR     = p.Results.TR;

if ~exist(OUTDIR,'dir'), mkdir(OUTDIR); end

if nargin < 1 || isempty(taskInput)
    taskFiles = defaultTaskFiles;

elseif isstruct(taskInput)
    taskFiles = taskInput;

else
    files = string(taskInput);
    taskFiles = struct();

    for i = 1:numel(files)
        f = files(i);

        % If not found, try inside INDIR
        if ~isfile(f)
            f2 = fullfile(INDIR, f);
            if isfile(f2), f = f2; end
        end

        [~, base, ext] = fileparts(f);
        assert(strcmpi(ext, '.mat'), 'Not a .mat file: %s', f);

        taskName = regexprep(base, '^observed_', '');
        taskName = upper(taskName);
        taskName = matlab.lang.makeValidName(taskName);

        taskFiles.(taskName) = f;
    end
end

%% ---- PROCESS ALL TASKS ----
taskNames = fieldnames(taskFiles);
CrossTask = [];

for ti = 1:numel(taskNames)
    task = taskNames{ti};
    fmat = taskFiles.(task);
    if ~isfile(fmat)
        warning('File missing for %s: %s (skipping)', task, fmat);
        continue;
    end
    fprintf('== %s ==\n', task);

    % Load required variables
    S = load(fmat);
    required = {'TA_mean','TA_se','PMZ_mean','PMZ_se','Haufe_mean','Haufe_std','W_obs','contrast_labels'};
    for q = 1:numel(required)
        assert(isfield(S, required{q}), 'Missing %s in %s', required{q}, fmat);
    end

    TA     = S.TA_mean;
    TA_se  = S.TA_se;
    PMZ    = S.PMZ_mean;
    PMZ_se = S.PMZ_se;

    Hmean = S.Haufe_mean;  % [R x T x nCon]
    Hstd  = S.Haufe_std;   % [R x T x nCon]

    % Calculate Weight Statistics from raw W_obs [Rg x T x nCon x Boot]
    Wmean = mean(S.W_obs, 4);
    Wstd  = std(S.W_obs, 0, 4);

    labels = S.contrast_labels(:);
    T = size(TA,1); nCon = size(TA,2);
    R = size(Hmean,1);

    tvec = (0:T-1) * TR;

    % ---- Stage 1: Temporal Summaries ----
    Sum = table();
    Sum.Contrast       = string(labels);
    Sum.PeakTA         = nan(nCon,1);
    Sum.PeakTA_t       = nan(nCon,1);
    Sum.PeakPMZ        = nan(nCon,1);
    Sum.PeakPMZ_t      = nan(nCon,1);
    Sum.PeakAbsPMZ     = nan(nCon,1);
    Sum.PeakAbsPMZ_t   = nan(nCon,1);

    % Store indices for spatial lookup later
    peak_indices = nan(nCon,1);

    for c = 1:nCon
        [pkTA, it]     = max(TA(:,c), [], 'omitnan');
        [pkPMZ, ipmz]  = max(PMZ(:,c), [], 'omitnan');
        [pkAPMZ, iap]  = max(abs(PMZ(:,c)), [], 'omitnan');

        if isempty(it)   || isnan(it),   it   = 1; end
        if isempty(ipmz) || isnan(ipmz), ipmz = 1; end
        if isempty(iap)  || isnan(iap),  iap  = 1; end

        Sum.PeakTA(c)       = pkTA;
        Sum.PeakTA_t(c)     = tvec(it);
        Sum.PeakPMZ(c)      = pkPMZ;
        Sum.PeakPMZ_t(c)    = tvec(ipmz);
        Sum.PeakAbsPMZ(c)   = pkAPMZ;
        Sum.PeakAbsPMZ_t(c) = tvec(iap);

        peak_indices(c) = ipmz; % Store PMZ peak index for spatial lookup
    end

    % ---- Stage 2: Spatial Analysis (Haufe vs Weights) ----

    % A. Calculate Z-scores (Mean / SD across bootstraps)
    H_Z = Hmean ./ (Hstd + eps);
    W_Z = Wmean ./ (Wstd + eps);

    % B. Identify Best Contrast (Highest PMZ Peak)
    [~, best_c_idx] = max(Sum.PeakPMZ);
    best_peak_t     = peak_indices(best_c_idx);
    best_con_label  = labels{best_c_idx};

    % C. Extract Spatial Stats for the Best Contrast at its Peak Time
    HZ_at_peak = H_Z(:, best_peak_t, best_c_idx);
    WZ_at_peak = W_Z(:, best_peak_t, best_c_idx);
    H_raw      = Hmean(:, best_peak_t, best_c_idx);

    % D. Sort by Absolute Haufe Z-score (Reliability)
    [~, sort_idx] = sort(abs(HZ_at_peak), 'descend');
    top50_idx = sort_idx(1:min(50, R));

    % E. Create Spatial Stats Table
    SpatialStats = table();
    SpatialStats.Rank = (1:numel(top50_idx))';
    SpatialStats.RegionID = top50_idx;
    SpatialStats.Haufe_Z = HZ_at_peak(top50_idx);
    SpatialStats.Weight_Z = WZ_at_peak(top50_idx);
    SpatialStats.Haufe_Raw = H_raw(top50_idx);

    % ---- Write CSVs ----
    taskSafe = lower(task);
    writetable(Sum, fullfile(OUTDIR, sprintf('%s_summary_contrasts.csv', taskSafe)));

    spatial_csv = fullfile(OUTDIR, sprintf('%s_spatial_stats_best_contrast.csv', taskSafe));
    writetable(SpatialStats, spatial_csv);
    fprintf('   Spatial stats for best contrast (%s) saved.\n', best_con_label);

    % ---- Plots: Timecourses with RIBBONS ----
    plot_timecourses(tvec, TA, TA_se, labels, sprintf('%s — TA', task), ...
        fullfile(OUTDIR, sprintf('%s_TA_timecourses.png', taskSafe)), 0.5, ...
        'Accuracy (± SE)', [0.3 1.0]);

    plot_timecourses(tvec, PMZ, PMZ_se, labels, sprintf('%s — PMZ', task), ...
        fullfile(OUTDIR, sprintf('%s_PMZ_timecourses.png', taskSafe)), 0, ...
        'Paired standardized margin (± SE)', []);

    % ---- Plots: Top 50 Haufe Heatmap (Time x Region) ----
    Top50_Timecourse = squeeze(H_Z(top50_idx, :, best_c_idx));

    plot_region_heatmap(Top50_Timecourse, tvec, top50_idx, ...
        sprintf('%s Top 50 Regions (Haufe Z) - %s', task, best_con_label), ...
        fullfile(OUTDIR, sprintf('%s_Haufe_Top50.png', taskSafe)));

    % ---- Cross-task summary pickers ----
    [bestPMZ, icPMZ] = max(Sum.PeakPMZ, [], 'omitnan');

    CTask = struct( ...
        'Task', task, ...
        'BestContrast', Sum.Contrast(icPMZ), ...
        'PeakPMZ', bestPMZ, ...
        'PeakTime', Sum.PeakPMZ_t(icPMZ), ...
        'TR', TR ...
    );
    CrossTask = [CrossTask; CTask]; %#ok<AGROW>
end

%% ---- Cross-task CSV summary ----
if ~isempty(CrossTask)
    Tcross = struct2table(CrossTask);
    writetable(Tcross, fullfile(OUTDIR, 'cross_task_best_contrasts.csv'));
end

fprintf('Done. Outputs in: %s\n', OUTDIR);
end

% ================== helpers ==================

function plot_timecourses(tvec, Y, E, labels, ttl, outpng, baseline, ylab, ylims)
figure('Color','w','Position',[100 100 1100 420]);
hold on
if ~isempty(baseline)
    yline(baseline, '--', sprintf('Reference = %.1f', baseline), 'Color',[.3 .3 .3], 'LineWidth', 1.2);
end

nCon = size(Y, 2);
colors = lines(nCon);

for c = 1:nCon
    curve = Y(:, c)';
    err   = E(:, c)';
    x_fill = [tvec, fliplr(tvec)];
    y_fill = [curve + err, fliplr(curve - err)];
    fill(x_fill, y_fill, colors(c,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(tvec, curve, 'Color', colors(c,:), 'LineWidth', 2, 'DisplayName', labels{c});
end

xlabel('Time (s)'); ylabel(ylab);
title(ttl, 'Interpreter','none');
if numel(labels) <= 20, legend('Interpreter','none','Location','eastoutside'); end
grid on; xlim([tvec(1) tvec(end)]);
if ~isempty(ylims), ylim(ylims); end
exportgraphics(gcf, outpng, 'Resolution', 200);
close;
end

function plot_region_heatmap(Data, tvec, regionIDs, ttl, outpng)
% Data: [Regions x Time]
figure('Color','w','Position',[100 100 1000 800]);

y_labels = arrayfun(@(x) sprintf('Reg %d', x), regionIDs, 'UniformOutput', false);

imagesc(tvec, 1:length(regionIDs), Data);
xlabel('Time (s)'); ylabel('Ranked Region ID');
title(ttl, 'Interpreter','none');
cb = colorbar; cb.Label.String = 'Haufe Z-Score';

colormap(parula);
clim([-3 3]);

yticks(1:length(regionIDs));
yticklabels(y_labels);
set(gca, 'FontSize', 8);

exportgraphics(gcf, outpng, 'Resolution', 200);
close;
end
