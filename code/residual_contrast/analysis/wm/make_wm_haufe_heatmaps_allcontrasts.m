%% make_wm_haufe_heatmaps_allcontrasts.m
clear; clc;

plot_labels = { ...
    'Body vs Face', ...
    'Body vs Place', ...
    'Body vs Tool', ...
    'Face vs Place', ...
    'Face vs Tool', ...
    'Place vs Tool'};

obsDir = fullfile(pwd, 'data', 'observed_residual_results');

files = { ...
    'observed_wm_interaction_chrf.mat', ...
    'observed_wm_interaction_chrfderiv.mat', ...
    'observed_wm_interaction_shrf.mat'};

model_names = {'cHRF', 'cHRF+deriv', 'sHRF'};
TR = 0.72;

scaleMode = 'all';   % 'all' or 'contrast'
clipPct = 99;

outFile = fullfile(pwd, 'figures', 'WM_Haufe_heatmaps_allcontrasts.png');
if ~exist(fileparts(outFile), 'dir')
    mkdir(fileparts(outFile));
end

for m = 1:3
    D = load(fullfile(obsDir, files{m}), 'Haufe_mean', 'contrast_labels');
    M(m).name = model_names{m};
    M(m).H = D.Haufe_mean;   % [489 x 41 x 6]
    M(m).contrast_labels = plot_labels;
end

contrast_labels = M(1).contrast_labels;
[R, T, C] = size(M(1).H);
tsec = (0:T-1) * TR;

switch lower(scaleMode)
    case 'all'
        allvals = [];
        for m = 1:3
            tmp = M(m).H;
            allvals = [allvals; tmp(:)];
        end
        mx = prctile(abs(allvals), clipPct);
        CLIM_ALL = [-mx mx];

    case 'contrast'
        CLIM_CONTRAST = zeros(C, 2);
        for c = 1:C
            vals = [];
            for m = 1:3
                tmp = M(m).H(:,:,c);
                vals = [vals; tmp(:)];
            end
            mx = prctile(abs(vals), clipPct);
            CLIM_CONTRAST(c,:) = [-mx mx];
        end

    otherwise
        error('scaleMode must be ''all'' or ''contrast''.');
end

f = figure('Color', 'w', 'Position', [100 50 1800 2600]);
tl = tiledlayout(C, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

for c = 1:C
    for m = 1:3
        nexttile;
        imagesc(tsec, 1:R, M(m).H(:,:,c));
        axis xy;
        set(gca, 'FontSize', 13);

        if strcmpi(scaleMode, 'all')
            clim(CLIM_ALL);
        else
            clim(CLIM_CONTRAST(c,:));
        end

        colormap(gca, bluewhitered_local(256));

        % Show seconds only on bottom row
        if c == C
            xlabel('Seconds', 'FontSize', 13);
        else
            set(gca, 'XTickLabel', []);
        end

        % Larger model titles
        if c == 1
            title(M(m).name, 'FontWeight', 'bold', 'FontSize', 16);
        end

        if m == 1
            ylabel({strrep(contrast_labels{c}, '_', '\_'), 'Region'}, ...
                'Interpreter', 'tex', 'FontSize', 14);
        else
            ylabel('Region', 'FontSize', 14);
        end
    end
end

cb = colorbar;
cb.Layout.Tile = 'east';
cb.FontSize = 14;
cb.Label.String = 'Haufe weight';
cb.Label.FontSize = 14;

exportgraphics(f, outFile, 'Resolution', 300);
fprintf('Saved: %s\n', outFile);

function cmap = bluewhitered_local(n)
    if nargin < 1
        n = 256;
    end

    n1 = floor(n/2);
    n2 = n - n1;

    blue = [0 0 1];
    white = [1 1 1];
    red = [1 0 0];

    cmap1 = [linspace(blue(1),white(1),n1)' ...
             linspace(blue(2),white(2),n1)' ...
             linspace(blue(3),white(3),n1)'];

    cmap2 = [linspace(white(1),red(1),n2)' ...
             linspace(white(2),red(2),n2)' ...
             linspace(white(3),red(3),n2)'];

    cmap = [cmap1; cmap2];
end