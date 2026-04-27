baseDir = 'C:\Users\Ryan_\Documents\my_Git\fMRI\data\observed_residual_results';
outdir  = fullfile(baseDir, 'language_simple_analysis');
if ~exist(outdir, 'dir'), mkdir(outdir); end

files = { ...
    fullfile(baseDir,'observed_language_simple_chrf.mat'), ...
    fullfile(baseDir,'observed_language_simple_chrfderiv.mat'), ...
    fullfile(baseDir,'observed_language_simple_shrf.mat')};

modelNames = {'cHRF','cHRF+deriv','sHRF'};
%% PairDiff_mean
data = cell(1,3);
for m = 1:3
    S = load(files{m});
    data{m}.PMZ_mean = S.PMZ_mean;
end

time_s = (0:size(data{1}.PMZ_mean,1)-1) * 0.72;

titles = { ...
    'Present Phase: Story vs Math', ...
    'Question Phase: Story vs Math', ...
    'Response Phase: Story vs Math'};

fig = figure('Color','w','Position',[100 100 1350 300]);
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

for k = 1:3
    nexttile; hold on

    for m = 1:3
        plot(time_s, data{m}.PMZ_mean(:,k), 'LineWidth', 1.8);
    end

    yline(0,'k--','LineWidth',1);

    title(titles{k}, 'FontSize', 10, 'FontWeight', 'bold');
    xlabel('Seconds');
    ylabel('Paired Score Difference');
    xlim([0 29]);
    grid off
    box off
    set(gca, 'FontSize', 9, 'LineWidth', 0.8);
end

legend(modelNames, 'Location','northwest', 'Box','off', 'FontSize',8);

exportgraphics(fig, fullfile(outdir, 'language_story_math_pairdiff_3panels.png'), 'Resolution', 300);

%% TA
data = cell(1,3);
for m = 1:3
    S = load(files{m});
    data{m}.TA_mean = S.TA_mean;
end

time_s = (0:size(data{1}.TA_mean,1)-1) * 0.72;

titles = { ...
    'Present Phase: Story vs Math', ...
    'Question Phase: Story vs Math', ...
    'Response Phase: Story vs Math'};

fig = figure('Color','w','Position',[100 100 1350 300]);
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

for k = 1:3
    nexttile; hold on

    for m = 1:3
        plot(time_s, data{m}.TA_mean(:,k), 'LineWidth', 1.8);
    end

    yline(0.5,'k--','LineWidth',1);

    title(titles{k}, 'FontSize', 10, 'FontWeight', 'bold');
    xlabel('Seconds');
    ylabel('Temporal Accuracy');
    xlim([0 29]);
    grid off
    box off
    set(gca, 'FontSize', 9, 'LineWidth', 0.8);
end

legend(modelNames, 'Location','northwest', 'Box','off', 'FontSize',8);

exportgraphics(fig, fullfile(outdir, 'language_story_math_3panels.png'), 'Resolution', 300);

%% LANGUAGE: mean |Haufe| overview for 3 between-task contrasts only
% M(m).H is [regions x time x contrasts]
% m = 1:cHRF, 2:cHRFderiv, 3:sHRF

nModels = numel(M);
model_names = {M.name};
all_labels = M(1).contrast_labels;

% colors
if ~exist('cols','var') || size(cols,1) < nModels
    cols = lines(nModels);
end

% fallback t-axis if tsec is not already defined
[~, T, ~] = size(M(1).H);
if ~exist('tsec','var') || numel(tsec) ~= T
    warning('tsec not found or wrong length; using 1:T instead.');
    tsec = 1:T;
end

%% Find the 3 desired between-task contrasts
want_pairs = {
    'StoryPresent',  'MathPresent';
    'StoryQuestion', 'MathQuestion';
    'StoryResponse', 'MathResponse'
};

plot_titles = {
    'Present phase: Story vs Math'
    'Question phase: Story vs Math'
    'Response phase: Story vs Math'
};

keep = zeros(1,3);

for k = 1:3
    a = want_pairs{k,1};
    b = want_pairs{k,2};

    idx = find( ...
        cellfun(@(x) contains(x,a,'IgnoreCase',true) && contains(x,b,'IgnoreCase',true), all_labels) ...
    );

    if isempty(idx)
        error('Could not find contrast containing both "%s" and "%s".', a, b);
    elseif numel(idx) > 1
        warning('Multiple matches found for %s vs %s; using the first one.', a, b);
        idx = idx(1);
    end

    keep(k) = idx;
end

% subset contrasts
C = numel(keep);
contrast_labels = all_labels(keep);

%% Mean absolute Haufe over time
meanAbsHaufe = cell(1, nModels);
for m = 1:nModels
    % subset to selected contrasts only
    Hsub = M(m).H(:,:,keep);                     % [R x T x 3]
    meanAbsHaufe{m} = squeeze(mean(abs(Hsub), 1, 'omitnan'));   % [T x 3]
end

%% Plot overview
f3 = figure('Color','w','Position',[100 100 1350 420]);
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

for c = 1:C
    nexttile; hold on;

    for m = 1:nModels
        plot(tsec, meanAbsHaufe{m}(:,c), 'LineWidth', 2, 'Color', cols(m,:));
    end

    xline(0,'k:');   % keep onset line only

    xlabel('Seconds');
    ylabel('Mean |Haufe|');
    title(plot_titles{c});

    if c == 1
        legend(model_names, 'Location','best', 'Box','off');
    end
end

exportgraphics(f3, 'LANGUAGE_HaufeMagnitude_betweenTaskOnly.png', 'Resolution', 300);