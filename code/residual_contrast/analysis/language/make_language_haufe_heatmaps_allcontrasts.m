clear; clc;

plot_labels = { ...
    'Story vs Math: Stimulus', ...
    'Story vs Math: Question', ...
    'Story vs Math: Response'};

phase_tags = {'stimulus', 'question', 'response'};

obsDir = fullfile(pwd, 'data', 'observed_residual_results');

files = { ...
    'observed_language_simple_chrf.mat', ...
    'observed_language_simple_chrfderiv.mat', ...
    'observed_language_simple_shrf.mat'};

model_names = {'cHRF', 'cHRF+deriv', 'sHRF'};
TR = 0.72;

scaleMode = 'all';   % 'all' or 'contrast'
clipPct = 99;

outDir = fullfile(pwd, 'figures');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% Load data
for m = 1:3
    D = load(fullfile(obsDir, files{m}), 'Haufe_mean', 'contrast_labels');

    M(m).name = model_names{m};
    M(m).H = D.Haufe_mean;   % expected: [R x T x 3]

    if isfield(D, 'contrast_labels') && numel(D.contrast_labels) == size(D.Haufe_mean, 3)
        M(m).contrast_labels = D.contrast_labels;
    else
        M(m).contrast_labels = plot_labels;
    end
end

contrast_labels = M(1).contrast_labels;
[R, T, C] = size(M(1).H);
tsec = (0:T-1) * TR;

%% Scaling
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

%% Make one figure per phase
for c = 1:C

    f = figure('Color', 'w', 'Position', [100 200 1800 400]);
    tl = tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

    for m = 1:3
        ax = nexttile;
        imagesc(tsec, 1:R, M(m).H(:,:,c));
        axis xy;
        pbaspect([4 1 1]);   % make heatmap itself shorter
        set(gca, 'FontSize', 13);

        if strcmpi(scaleMode, 'all')
            clim(CLIM_ALL);
        else
            clim(CLIM_CONTRAST(c,:));
        end

        colormap(gca, bluewhitered_local(256));

        xlabel('Seconds', 'FontSize', 13);

        if m == 1
            ylabel('Region', 'FontSize', 14);
        else
            ylabel('');
        end

        % Model labels above each heatmap
        text(ax, 0.5, 1.03, M(m).name, ...
            'Units', 'normalized', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontWeight', 'bold', ...
            'FontSize', 16);
    end

    cb = colorbar;
    cb.Layout.Tile = 'east';
    cb.FontSize = 14;
    cb.Label.String = 'Haufe weight';
    cb.Label.FontSize = 14;

    outFile = fullfile(outDir, ...
        sprintf('LANGUAGE_Haufe_heatmap_%s.png', phase_tags{c}));

    exportgraphics(f, outFile, 'Resolution', 300);
    fprintf('Saved: %s\n', outFile);
end

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