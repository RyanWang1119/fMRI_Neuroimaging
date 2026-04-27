%%
plot_labels = { ...
    'Body vs Face', ...
    'Body vs Place', ...
    'Body vs Tool', ...
    'Face vs Place', ...
    'Face vs Tool', ...
    'Place vs Tool'};

files = { ...
    'observed_wm_interaction_chrf.mat', ...
    'observed_wm_interaction_chrfderiv.mat', ...
    'observed_wm_interaction_shrf.mat'};

model_names = {'cHRF', 'cHRF+deriv', 'sHRF'};
TR = 0.72;

for m = 1:3
    D = load(files{m});

    M(m).name            = model_names{m};
    M(m).TA_mean         = D.TA_mean;          % [41 x 6]
    M(m).PMZ_mean        = D.PMZ_mean;         % [41 x 6]
    M(m).Haufe_mean      = D.Haufe_mean;       % [489 x 41 x 6]
    M(m).Haufe_obs       = D.Haufe_obs;        % [489 x 41 x 6 x 100]
    M(m).contrast_labels = plot_labels;
end

T = size(M(1).TA_mean, 1);
C = size(M(1).TA_mean, 2);
Rg = size(M(1).Haufe_mean, 1);

tsec = (0:T-1) * TR;
contrast_labels = M(1).contrast_labels;

cols = lines(3);

f1 = figure('Color','w','Position',[100 100 1300 750]);
tl = tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

for c = 1:C
    nexttile; hold on;

    for m = 1:3
        plot(tsec, M(m).TA_mean(:,c), 'LineWidth', 2, 'Color', cols(m,:));
    end

    xline(0,'k:');
    xlabel('Seconds');
    ylabel('Temporal Accuracy');
    title(strrep(contrast_labels{c}, '_', '\_'));

    if c == 1
        legend(model_names, 'Location','best', 'Box','off');
    end
end

exportgraphics(f1, 'WM_TA_overview.png', 'Resolution', 300);

%%
f2 = figure('Color','w','Position',[100 100 1300 750]);
tl = tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

for c = 1:C
    nexttile; hold on;

    for m = 1:3
        plot(tsec, M(m).PMZ_mean(:,c), 'LineWidth', 2, 'Color', cols(m,:));
    end

    xline(0,'k:');
    xlabel('Seconds');
    ylabel('Paired score difference');
    title(strrep(contrast_labels{c}, '_', '\_'));

    if c == 1
        legend(model_names, 'Location','best', 'Box','off');
    end
end

exportgraphics(f2, 'WM_PMZ_overview.png', 'Resolution', 300);

%%
refModel = 1;   % cHRF
halfwin  = 1;   % peak +/- 1 TR

peak_idx = zeros(C,1);
win = cell(C,1);

for c = 1:C
    [~, peak_idx(c)] = max(M(refModel).TA_mean(:,c));
    win{c} = max(1, peak_idx(c)-halfwin) : min(T, peak_idx(c)+halfwin);
end

Zavg = cell(C,3);
allvals = [];

for c = 1:C
    for m = 1:3
        z = mean(M(m).Haufe_mean(:, win{c}, c), 2, 'omitnan');   % [489 x 1]
        Zavg{c,m} = z;
        allvals = [allvals; z(:)];
    end
end

% common symmetric color scale across all 18 maps
mx = prctile(abs(allvals), 99);
clim = [-mx mx];

for m = 1:3
    % Haufe_mean: [489 x 41 x 6]
    meanAbsHaufe{m} = squeeze(mean(abs(M(m).Haufe_mean), 1, 'omitnan'));  % [41 x 6]
end

f3 = figure('Color','w','Position',[100 100 1300 750]);
tl = tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

for c = 1:C
    nexttile; hold on;

    for m = 1:3
        plot(tsec, meanAbsHaufe{m}(:,c), 'LineWidth', 2, 'Color', cols(m,:));
    end

    xline(0,'k:');
    xlabel('Seconds');
    ylabel('Mean |Haufe|');
    title(strrep(contrast_labels{c}, '_', '\_'));

    if c == 1
        legend(model_names, 'Location','best', 'Box','off');
    end
end

exportgraphics(f3, 'WM_HaufeMagnitude_overview.png', 'Resolution', 300);

%%
for m = 1:3
    H = M(m).Haufe_obs;   % [489 x 41 x 6 x 100]

    signCons = zeros(T, C);

    for c = 1:C
        for t = 1:T
            X = squeeze(H(:,t,c,:));   % [489 x 100]

            % sign consistency per region: abs(mean(sign over boots)))
            sc_region = abs(mean(sign(X), 2, 'omitnan'));   % [489 x 1], range 0..1

            signCons(t,c) = mean(sc_region, 'omitnan');
        end
    end

    M(m).signCons = signCons;
end

f4 = figure('Color','w','Position',[100 100 1300 750]);
tl = tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

for c = 1:C
    nexttile; hold on;

    for m = 1:3
        plot(tsec, M(m).signCons(:,c), 'LineWidth', 2, 'Color', cols(m,:));
    end

    xline(0,'k:');
    ylim([0 1]);
    xlabel('Seconds');
    ylabel('Mean sign consistency');
    title(strrep(contrast_labels{c}, '_', '\_'));

    if c == 1
        legend(model_names, 'Location','best', 'Box','off');
    end
end

exportgraphics(f4, 'WM_HaufeStability_overview.png', 'Resolution', 300);

