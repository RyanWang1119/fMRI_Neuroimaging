function visualize_spatial_patterns(taskName, targetContrast, targetTime, varargin)
% visualize_spatial_patterns('wm_interaction_chrf')
% visualize_spatial_patterns('wm_interaction_chrf','contrast_idx',1,'time_window',8:10)
%
% Supports:
%   - auto contrast/time selection by PMZ
%   - specific contrast by name or index
%   - specific time point or averaged time window
%   - export to file
%   - fixed color limits for cross-model comparisons
%
% For time-window plots, the function averages Haufe_obs across the selected
% window within each bootstrap, then computes mean and SE across bootstraps.
% This is preferable to averaging precomputed Z maps across time.

if nargin < 2, targetContrast = []; end
if nargin < 3, targetTime = []; end

p = inputParser;
addParameter(p, 'contrast_idx', [], @(x)isnumeric(x) && isscalar(x));
addParameter(p, 'time_window', [], @isnumeric);
addParameter(p, 'outfile', '', @(x)ischar(x) || isstring(x));
addParameter(p, 'clim', [], @isnumeric);
addParameter(p, 'title_str', '', @(x)ischar(x) || isstring(x));
parse(p, varargin{:});

contrast_idx_user = p.Results.contrast_idx;
time_window       = p.Results.time_window;
outfile           = char(p.Results.outfile);
clim              = p.Results.clim;
title_str_user    = char(p.Results.title_str);

% 1) Locate file and load needed fields
obsDir  = fullfile(pwd, 'data', 'observed_residual_results');
oldDir  = fullfile(pwd, 'data', 'task_residual');
obsFile = fullfile(obsDir, ['observed_' lower(taskName) '.mat']);
oldFile = fullfile(oldDir, ['result_' lower(taskName) '.mat']);

if isfile(obsFile)
    fmat = obsFile;
elseif isfile(oldFile)
    fmat = oldFile;
else
    error('File not found. Tried: %s and %s', obsFile, oldFile);
end

fprintf('Loading %s...\n', fmat);

% Load Haufe_obs if available, because it is best for window averaging
S = load(fmat, 'Haufe_mean','Haufe_std','Haufe_obs','PMZ_mean','contrast_labels','nBoot');

% 2) Pick contrast
if ~isempty(contrast_idx_user)
    con_idx = contrast_idx_user;
    if con_idx < 1 || con_idx > numel(S.contrast_labels)
        error('contrast_idx out of range.');
    end
elseif isempty(targetContrast)
    [max_pmz, linear_idx] = max(S.PMZ_mean(:));
    [peak_t_auto, con_idx] = ind2sub(size(S.PMZ_mean), linear_idx);
    fprintf('Auto-selected strongest contrast by PMZ: %s\n', S.contrast_labels{con_idx});
else
    con_idx = find(strcmp(S.contrast_labels, targetContrast), 1);
    if isempty(con_idx)
        fprintf('Exact match not found for "%s". Available contrasts:\n', targetContrast);
        disp(S.contrast_labels);
        error('Invalid contrast name.');
    end
    [max_pmz, peak_t_auto] = max(S.PMZ_mean(:, con_idx));
end

real_label = S.contrast_labels{con_idx};

% 3) Pick time or time window
if ~isempty(time_window)
    tw = unique(time_window(:))';
    fprintf('Using user-defined time window: %s\n', mat2str(tw));
    time_label = sprintf('t=%d-%d', tw(1), tw(end));
elseif isempty(targetTime)
    if ~exist('peak_t_auto','var')
        [max_pmz, peak_t_auto] = max(S.PMZ_mean(:, con_idx));
    end
    tw = peak_t_auto;
    fprintf('Auto-selected peak time: t=%d (PMZ: %.3f)\n', peak_t_auto, max_pmz);
    time_label = sprintf('t=%d', peak_t_auto);
else
    tw = targetTime;
    fprintf('Using user-defined time: t=%d (PMZ: %.3f)\n', targetTime, S.PMZ_mean(targetTime, con_idx));
    time_label = sprintf('t=%d', targetTime);
end

% 4) Build Haufe Z-like map
% Preferred route: use Haufe_obs, average within bootstrap across the window,
% then compute mean and SE across bootstraps.
if isfield(S, 'Haufe_obs') && ~isempty(S.Haufe_obs)
    % Haufe_obs: [region x time x contrast x boot]
    Hboot = squeeze(S.Haufe_obs(:, tw, con_idx, :));  % region x time/window x boot

    if isvector(Hboot)
        % edge case: single region/time/boot, force dims
        Hboot = reshape(Hboot, [numel(Hboot), 1, 1]);
    end

    if ndims(Hboot) == 2
        % region x boot (single time)
        Hwin_boot = Hboot;
    else
        % average over selected time window within each bootstrap
        Hwin_boot = squeeze(mean(Hboot, 2, 'omitnan'));   % region x boot
    end

    H_mean = mean(Hwin_boot, 2, 'omitnan');
    H_std  = std(Hwin_boot, 0, 2, 'omitnan');

    if isfield(S,'nBoot') && ~isempty(S.nBoot) && S.nBoot > 1
        H_se = H_std ./ sqrt(S.nBoot);
    else
        warning('nBoot missing; falling back to SD.');
        H_se = H_std;
    end

else
    % fallback route if Haufe_obs is not available
    H_mean = mean(S.Haufe_mean(:, tw, con_idx), 2, 'omitnan');
    H_std  = mean(S.Haufe_std(:, tw, con_idx), 2, 'omitnan');

    if isfield(S,'nBoot') && ~isempty(S.nBoot) && S.nBoot > 1
        H_se = H_std ./ sqrt(S.nBoot);
    else
        warning('nBoot missing; falling back to SD.');
        H_se = H_std;
    end
end

Zvals = H_mean ./ max(H_se, eps);

% 5) Map to brain using Canlab 2018 atlas
try
    atlas_obj = load_atlas('canlab2018');
catch
    error('Could not load atlas. Ensure CanLabCore is on your path.');
end

if numel(Zvals) ~= numel(atlas_obj.labels)
    warning('Region count mismatch: Data=%d, Atlas=%d.', numel(Zvals), numel(atlas_obj.labels));
end

r = atlas2region(atlas_obj);
for pidx = 1:numel(r)
    if pidx <= numel(Zvals)
        v = Zvals(pidx);
    else
        v = 0;
    end
    r(pidx).Z = v + zeros(size(r(pidx).Z));
end

% 6) Visualize
% Do NOT rely on a fresh figure handle here, because orthviews may draw
% into the SPM graphics window instead of the current figure.

orthviews(r);
drawnow;

% Try to get the actual figure used by orthviews / SPM
fig = [];
try
    fig = spm_figure('FindWin', 'Graphics');
catch
end

% Fallback if SPM figure handle is not available
if isempty(fig) || ~ishandle(fig)
    fig = gcf;
end

if ~isempty(clim)
    ax = findall(fig, 'Type', 'Axes');
    for k = 1:numel(ax)
        try
            caxis(ax(k), clim);
        end
    end
end

if isempty(title_str_user)
    title_str = sprintf('%s: %s (%s) — Haufe Z', taskName, real_label, time_label);
else
    title_str = title_str_user;
end

sgtitle(fig, strrep(title_str, '_', '\_'), 'Interpreter', 'none');
fprintf('Displaying: %s\n', title_str);

if ~isempty(outfile)
    drawnow;
    pause(0.2);  % helps orthviews finish rendering

    try
        exportgraphics(fig, outfile, 'Resolution', 300);
    catch
        % fallback: capture pixels from the displayed figure
        fr = getframe(fig);
        imwrite(fr.cdata, outfile);
    end
end