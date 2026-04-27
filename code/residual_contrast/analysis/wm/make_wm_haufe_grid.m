%% make_wm_haufe_grid.m
% Creates a 6 x 3 WM Haufe grid:
%   rows    = 6 WM contrasts
%   columns = cHRF, cHRF+deriv, sHRF
%
% Assumes your updated visualize_spatial_patterns.m is on the MATLAB path
% and supports:
%   'contrast_idx', 'time_window', 'outfile', 'clim', 'title_str'
%
% It uses one matched time window per contrast, chosen from the reference
% model, then applies that same window to all 3 models.

clear; clc;

%% settings
tags = { ...
    'wm_interaction_chrf', ...
    'wm_interaction_chrfderiv', ...
    'wm_interaction_shrf'};

model_names = {'cHRF', 'cHRF+deriv', 'sHRF'};

TR = 0.72;

% choose which model defines the matched time window
refModel = 1;          % 1=cHRF, 2=cHRF+deriv, 3=sHRF

% choose which metric defines the window center
% valid: 'TA' or 'PMZ'
windowMetric = 'TA';

% use peak +/- halfwin TRs
halfwin = 1;

% color scale percentile for robust common clim
clipPct = 99;

% output folders
panelDir = fullfile(pwd, 'figures', 'wm_haufe_panels');
if ~exist(panelDir, 'dir')
    mkdir(panelDir);
end

outFigure = fullfile(pwd, 'figures', 'WM_Haufe_grid.png');
if ~exist(fileparts(outFigure), 'dir')
    mkdir(fileparts(outFigure));
end

% data folder expected by visualize_spatial_patterns
obsDir = fullfile(pwd, 'data', 'observed_residual_results');

%% load model summaries
M = struct([]);

for m = 1:numel(tags)
    fmat = fullfile(obsDir, ['observed_' tags{m} '.mat']);
    if ~isfile(fmat)
        error('Could not find file: %s', fmat);
    end

    S = load(fmat, ...
        'TA_mean', 'PMZ_mean', 'Haufe_obs', 'Haufe_mean', 'Haufe_std', ...
        'contrast_labels', 'nBoot');

    M(m).tag             = tags{m};
    M(m).name            = model_names{m};
    M(m).TA_mean         = S.TA_mean;          % [41 x 6]
    M(m).PMZ_mean        = S.PMZ_mean;         % [41 x 6]
    M(m).Haufe_obs       = S.Haufe_obs;        % [489 x 41 x 6 x 100]
    M(m).Haufe_mean      = S.Haufe_mean;       % [489 x 41 x 6]
    M(m).Haufe_std       = S.Haufe_std;        % [489 x 41 x 6]
    M(m).contrast_labels = S.contrast_labels;
    M(m).nBoot           = S.nBoot;
end

contrast_labels = M(1).contrast_labels;
T = size(M(1).TA_mean, 1);
C = size(M(1).TA_mean, 2);
tsec = (0:T-1) * TR;

%% choose matched time window per contrast from reference model
peak_idx = zeros(C, 1);
win = cell(C, 1);

for c = 1:C
    switch upper(windowMetric)
        case 'TA'
            y = M(refModel).TA_mean(:, c);
        case 'PMZ'
            y = M(refModel).PMZ_mean(:, c);
        otherwise
            error('windowMetric must be ''TA'' or ''PMZ''.');
    end

    [~, peak_idx(c)] = max(y);
    win{c} = max(1, peak_idx(c)-halfwin) : min(T, peak_idx(c)+halfwin);
end

%% compute a common color scale across all 18 window-averaged Haufe Z maps
allZ = [];

for c = 1:C
    for m = 1:numel(M)
        zvals = compute_windowed_haufe_z(M(m), c, win{c});
        allZ = [allZ; zvals(:)]; %#ok<AGROW>
    end
end

mx = prctile(abs(allZ), clipPct);
if mx <= 0 || isnan(mx)
    mx = max(abs(allZ));
end
clim = [-mx, mx];

fprintf('Using common clim = [%.3f, %.3f]\n', clim(1), clim(2));

%% export 18 individual panel PNGs
panel_png = cell(C, numel(M));

for c = 1:C
    for m = 1:numel(M)
        outpng = fullfile(panelDir, sprintf('wm_haufe_c%d_m%d.png', c, m));

        ttl = sprintf('%s | %s | %.2f-%.2f s', ...
            contrast_labels{c}, M(m).name, ...
            tsec(win{c}(1)), tsec(win{c}(end)));

        fprintf('Exporting panel: contrast %d / model %d -> %s\n', c, m, outpng);

        visualize_spatial_patterns(M(m).tag, [], [], ...
            'contrast_idx', c, ...
            'time_window', win{c}, ...
            'outfile', outpng, ...
            'clim', clim, ...
            'title_str', ttl);

        panel_png{c,m} = outpng;
    end
end

%% stitch 18 PNGs into one 6 x 3 figure
f = figure('Color', 'w', 'Position', [100 100 1500 2200]);
tl = tiledlayout(C, numel(M), 'Padding', 'compact', 'TileSpacing', 'compact');

for c = 1:C
    for m = 1:numel(M)
        nexttile;
        I = imread(panel_png{c,m});
        imshow(I);
        axis off;

        if c == 1
            title(M(m).name, 'Interpreter', 'none', 'FontWeight', 'bold');
        end

        if m == 1
            text(-0.08, 0.5, strrep(contrast_labels{c}, '_', '\_'), ...
                'Units', 'normalized', ...
                'Rotation', 90, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontWeight', 'bold', ...
                'Interpreter', 'tex');
        end
    end
end

title(tl, sprintf('WM matched-window Haufe maps (%s-selected windows)', windowMetric), ...
    'FontWeight', 'bold');

exportgraphics(f, outFigure, 'Resolution', 300);
fprintf('Saved stitched figure to: %s\n', outFigure);

%% optional: print chosen windows
fprintf('\nChosen windows by contrast:\n');
for c = 1:C
    fprintf('  %d) %s : TR %s  (%.2f-%.2f s)\n', ...
        c, contrast_labels{c}, mat2str(win{c}), ...
        tsec(win{c}(1)), tsec(win{c}(end)));
end

%% local helper
function Zvals = compute_windowed_haufe_z(Mm, con_idx, tw)
% Computes the same window-averaged Haufe Z-like map used for common clim.
% Preferred route:
%   average Haufe_obs within each bootstrap across the selected window,
%   then compute mean / SE across bootstraps.

    if isfield(Mm, 'Haufe_obs') && ~isempty(Mm.Haufe_obs)
        % Haufe_obs: [region x time x contrast x boot]
        Hboot = squeeze(Mm.Haufe_obs(:, tw, con_idx, :));

        if isvector(Hboot)
            Hboot = reshape(Hboot, [numel(Hboot), 1, 1]);
        end

        if ndims(Hboot) == 2
            % region x boot
            Hwin_boot = Hboot;
        else
            % region x time x boot -> region x boot
            Hwin_boot = squeeze(mean(Hboot, 2, 'omitnan'));
        end

        H_mean = mean(Hwin_boot, 2, 'omitnan');
        H_std  = std(Hwin_boot, 0, 2, 'omitnan');

        if ~isempty(Mm.nBoot) && Mm.nBoot > 1
            H_se = H_std ./ sqrt(Mm.nBoot);
        else
            H_se = H_std;
        end

    else
        % fallback if Haufe_obs is unavailable
        H_mean = mean(Mm.Haufe_mean(:, tw, con_idx), 2, 'omitnan');
        H_std  = mean(Mm.Haufe_std(:, tw, con_idx), 2, 'omitnan');

        if ~isempty(Mm.nBoot) && Mm.nBoot > 1
            H_se = H_std ./ sqrt(Mm.nBoot);
        else
            H_se = H_std;
        end
    end

    Zvals = H_mean ./ max(H_se, eps);
end