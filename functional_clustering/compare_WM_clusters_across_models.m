%% compare_WM_clusters_across_models.m
% Requires:
%   results_all   % forced-k=2 WM results
%   atlas_obj     % for parcel labels if you want recurring parcel tables
%
% Assumes:
%   results_all.sHRF{1:4}, results_all.cHRFderiv{1:4}, results_all.cHRF{1:4}
%   each result has:
%       .cluster_labels
%       .cluster_counts
%       .cluster_means
%       .keep_idx

clearvars -except results_all atlas_obj

model_names = {'sHRF','cHRFderiv','cHRF'};
contrast_names = {'body','faces','places','tools'};
even_thresh = 0.25;   % minority proportion threshold for "even split"

%% -----------------------------
%% 1) Minority-cluster characteristics
%% -----------------------------
rows = {};

for c = 1:numel(contrast_names)
    for m = 1:numel(model_names)
        res = results_all.(model_names{m}){c};

        assert(numel(res.cluster_counts) == 2, ...
            'This script expects forced-k=2 results.');

        [~, ord] = sort(res.cluster_counts(:), 'ascend'); % minority first
        kmin = ord(1);

        y = res.cluster_means(kmin, :);
        n1 = res.cluster_counts(kmin);
        nTot = sum(res.cluster_counts);

        [peak_val, peak_t] = max(y);
        [trough_val, trough_t] = min(y);

        auc_signed = sum(y);
        auc_abs = sum(abs(y));
        rms_val = sqrt(mean(y.^2));

        rows(end+1,:) = { ...
            model_names{m}, ...
            contrast_names{c}, ...
            kmin, ...
            n1, ...
            n1 / nTot, ...
            peak_t, ...
            peak_val, ...
            trough_t, ...
            trough_val, ...
            auc_signed, ...
            auc_abs, ...
            rms_val ...
            };
    end
end

Tminor_char = cell2table(rows, 'VariableNames', ...
    {'Model','Contrast','MinorityCluster','Size','Proportion', ...
     'PeakTime','PeakValue','TroughTime','TroughValue', ...
     'AUCSigned','AUCAbs','RMS'});

disp(Tminor_char)

%% -----------------------------
%% 2) Pairwise similarity of minority-cluster trajectories across models
%% -----------------------------
rows = {};

for c = 1:numel(contrast_names)
    for a = 1:numel(model_names)-1
        for b = a+1:numel(model_names)

            resA = results_all.(model_names{a}){c};
            resB = results_all.(model_names{b}){c};

            [~, ordA] = sort(resA.cluster_counts(:), 'ascend');
            [~, ordB] = sort(resB.cluster_counts(:), 'ascend');

            yA = resA.cluster_means(ordA(1), :);  % minority mean trajectory
            yB = resB.cluster_means(ordB(1), :);

            r = corr(yA(:), yB(:), 'rows', 'complete');
            rmse = sqrt(mean((yA(:) - yB(:)).^2));

            rows(end+1,:) = { ...
                contrast_names{c}, ...
                model_names{a}, ...
                model_names{b}, ...
                r, ...
                rmse ...
                };
        end
    end
end

Tminor_sim = cell2table(rows, 'VariableNames', ...
    {'Contrast','Model1','Model2','MinorityTrajectoryCorr','MinorityTrajectoryRMSE'});

disp(Tminor_sim)

%% -----------------------------
%% 3) Even-split recurrence analysis
%%    Only compare model pairs where BOTH have minority proportion >= even_thresh
%%    Use optimal 2-cluster matching by Jaccard overlap
%% -----------------------------
rows = {};
recurring_tables = struct();

for c = 1:numel(contrast_names)
    contrast_name = contrast_names{c};

    % determine which models qualify as "even split"
    qualifies = false(1, numel(model_names));
    props = nan(1, numel(model_names));

    for m = 1:numel(model_names)
        res = results_all.(model_names{m}){c};
        props(m) = min(res.cluster_counts) / sum(res.cluster_counts);
        qualifies(m) = props(m) >= even_thresh;
    end

    % only compare pairs where both qualify
    for a = 1:numel(model_names)-1
        for b = a+1:numel(model_names)
            if ~(qualifies(a) && qualifies(b))
                continue
            end

            resA = results_all.(model_names{a}){c};
            resB = results_all.(model_names{b}){c};

            match = match_two_cluster_solutions(resA, resB);

            % matched cluster overlaps
            for k = 1:2
                setA = match.setA{k};
                setB = match.setB_matched{k};

                inter = intersect(setA, setB);
                jac = safe_jaccard(setA, setB);

                yA = resA.cluster_means(k, :);
                yB = resB.cluster_means(match.permB(k), :);
                traj_corr = corr(yA(:), yB(:), 'rows', 'complete');

                rows(end+1,:) = { ...
                    contrast_name, ...
                    model_names{a}, ...
                    model_names{b}, ...
                    k, ...
                    numel(setA), ...
                    numel(setB), ...
                    numel(inter), ...
                    jac, ...
                    traj_corr ...
                    };

                % save recurring parcels for this matched cluster
                key = matlab.lang.makeValidName(sprintf('%s_%s_vs_%s_cluster%d', ...
                    contrast_name, model_names{a}, model_names{b}, k));

                recurring_tables.(key) = table(inter(:), atlas_obj.labels(inter(:)), ...
                    'VariableNames', {'ParcelIndex','ParcelName'});
            end
        end
    end
end

Teven_overlap = cell2table(rows, 'VariableNames', ...
    {'Contrast','Model1','Model2','MatchedCluster', ...
     'NModel1','NModel2','Intersection','Jaccard','TrajectoryCorr'});

disp(Teven_overlap)

%% -----------------------------
%% 4) Figures
%% -----------------------------

% Figure 1: minority-cluster residual trajectories across models
fig1 = figure('Color','w', 'Name', 'Minority cluster residual trajectories');
tiledlayout(2,2, 'Padding','compact', 'TileSpacing','compact');

for c = 1:numel(contrast_names)
    nexttile;
    hold on;

    for m = 1:numel(model_names)
        res = results_all.(model_names{m}){c};
        [~, ord] = sort(res.cluster_counts(:), 'ascend');
        kmin = ord(1);
        plot(res.cluster_means(kmin, :), 'LineWidth', 2, ...
            'DisplayName', sprintf('%s (n=%d)', model_names{m}, res.cluster_counts(kmin)));
    end

    yline(0, '--', 'LineWidth', 0.8);
    grid on;
    title(sprintf('%s', contrast_names{c}));
    xlabel('Time');
    ylabel('Mean residual');

    if c == 1
        legend('Location','best');
    end
end

sgtitle('Minority-cluster residual trajectories across HRF models');

% Figure 2: cluster-size proportions across models
fig2 = figure('Color','w', 'Name', 'Cluster size proportions');
tiledlayout(2,2, 'Padding','compact', 'TileSpacing','compact');

for c = 1:numel(contrast_names)
    nexttile;
    P = nan(numel(model_names), 2);

    for m = 1:numel(model_names)
        res = results_all.(model_names{m}){c};
        counts = sort(res.cluster_counts(:), 'descend');   % dominant, minority
        P(m,:) = counts ./ sum(counts);
    end

    bar(P, 'stacked', 'LineWidth', 1);
    xticks(1:numel(model_names));
    xticklabels(model_names);
    ylabel('Proportion');
    ylim([0 1]);
    title(sprintf('%s', contrast_names{c}));
    grid on;

    if c == 1
        legend({'Dominant','Minority'}, 'Location','southoutside');
    end
end

sgtitle('Cluster size proportions across HRF models');

% Figure 3: matched-cluster recurrence for even-split contrasts
if ~isempty(Teven_overlap)
    fig3 = figure('Color','w', 'Name', 'Matched cluster recurrence');
    tiledlayout(1,2, 'Padding','compact', 'TileSpacing','compact');

    nexttile;
    gscatter(categorical(Teven_overlap.Contrast), Teven_overlap.Jaccard, ...
        Teven_overlap.Model1 + " vs " + Teven_overlap.Model2, [], 'o', 8);
    ylabel('Jaccard overlap');
    title('Matched-cluster regional recurrence');
    grid on;

    nexttile;
    gscatter(categorical(Teven_overlap.Contrast), Teven_overlap.TrajectoryCorr, ...
        Teven_overlap.Model1 + " vs " + Teven_overlap.Model2, [], 'o', 8);
    ylabel('Trajectory correlation');
    title('Matched-cluster residual similarity');
    grid on;

    sgtitle(sprintf('Even-split contrasts (minority proportion >= %.2f)', even_thresh));
end

%% -----------------------------
%% Save outputs
%% -----------------------------
save('WM_cluster_crossmodel_analysis.mat', ...
    'Tminor_char', 'Tminor_sim', 'Teven_overlap', 'recurring_tables', ...
    'even_thresh', '-v7.3');

writetable(Tminor_char, 'WM_minority_cluster_characteristics.csv');
writetable(Tminor_sim, 'WM_minority_cluster_similarity.csv');

if ~isempty(Teven_overlap)
    writetable(Teven_overlap, 'WM_even_split_cluster_recurrence.csv');
end

%% -----------------------------
%% Helper functions
%% -----------------------------
function match = match_two_cluster_solutions(resA, resB)
    % Original parcel sets for each cluster
    setA = cell(2,1);
    setB = cell(2,1);

    for k = 1:2
        setA{k} = resA.keep_idx(resA.cluster_labels == k);
        setB{k} = resB.keep_idx(resB.cluster_labels == k);
    end

    % Jaccard matrix
    J = nan(2,2);
    for i = 1:2
        for j = 1:2
            J(i,j) = safe_jaccard(setA{i}, setB{j});
        end
    end

    % choose permutation of B maximizing total overlap
    score_direct = J(1,1) + J(2,2);
    score_swap   = J(1,2) + J(2,1);

    if score_direct >= score_swap
        permB = [1 2];
    else
        permB = [2 1];
    end

    setB_matched = cell(2,1);
    for k = 1:2
        setB_matched{k} = setB{permB(k)};
    end

    match.J = J;
    match.permB = permB;
    match.setA = setA;
    match.setB = setB;
    match.setB_matched = setB_matched;
end

function j = safe_jaccard(A, B)
    A = unique(A(:));
    B = unique(B(:));
    inter = numel(intersect(A,B));
    union = numel(union(A,B));
    if union == 0
        j = NaN;
    else
        j = inter / union;
    end
end