function visualize_spatial_patterns(taskName, targetContrast, targetTime)
% visualize_spatial_patterns('GAMBLING','loss_vs_win',10)
% Uses Haufe_SE (std/sqrt(nBoot)) to form Z-like stability maps.
% Auto-selection now uses PMZ (paired standardized margin), not TAFC.

if nargin < 2, targetContrast = []; end
if nargin < 3, targetTime = []; end

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
S = load(fmat, 'Haufe_mean','Haufe_std','PMZ_mean','contrast_labels','nBoot');

% 2) Pick contrast
if isempty(targetContrast)
    [max_pmz, linear_idx] = max(S.PMZ_mean(:));
    [peak_t_auto, con_idx] = ind2sub(size(S.PMZ_mean), linear_idx);
    fprintf('Auto-selected strongest contrast by PMZ: %s\n', S.contrast_labels{con_idx});
else
    con_idx = find(strcmp(S.contrast_labels, targetContrast));
    if isempty(con_idx)
        fprintf('Exact match not found for "%s". Available contrasts:\n', targetContrast);
        disp(S.contrast_labels);
        error('Invalid contrast name.');
    end
    [max_pmz, peak_t_auto] = max(S.PMZ_mean(:, con_idx));
end

% 3) Pick time
if isempty(targetTime)
    t = peak_t_auto;
    fprintf('Auto-selected peak time: t=%d (PMZ: %.3f)\n', t, max_pmz);
else
    t = targetTime;
    fprintf('Using user-defined time: t=%d (PMZ: %.3f)\n', t, S.PMZ_mean(t, con_idx));
end
real_label = S.contrast_labels{con_idx};

% 4) Build SE-based Z map
H_mean = S.Haufe_mean(:, t, con_idx);
H_std  = S.Haufe_std(:, t, con_idx);

if isfield(S,'nBoot') && ~isempty(S.nBoot) && S.nBoot > 1
    H_se = H_std ./ sqrt(S.nBoot);
else
    warning('nBoot missing; falling back to SD (conservative).');
    H_se = H_std;
end

Zvals = H_mean ./ max(H_se, eps);   % Haufe_Z (t-like stability)

% 5) Map to brain (Canlab 2018 atlas)
try
    atlas_obj = load_atlas('canlab2018');
catch
    error('Could not load atlas. Ensure CanLabCore is on your path.');
end

if numel(Zvals) ~= numel(atlas_obj.labels)
    warning('Region count mismatch: Data=%d, Atlas=%d.', numel(Zvals), numel(atlas_obj.labels));
end

r = atlas2region(atlas_obj);
for p = 1:numel(r)
    v = (p <= numel(Zvals)) * Zvals(min(p,numel(Zvals)));
    r(p).Z = v + zeros(size(r(p).Z));
end

% 6) Visualize
orthviews(r);
title_str = sprintf('%s: %s (t=%d) — Haufe Z [selected by PMZ]', taskName, real_label, t);
title(title_str);
fprintf('Displaying: %s\n', title_str);
end
