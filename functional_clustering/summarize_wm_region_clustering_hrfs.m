function T = summarize_wm_region_clustering_hrfs(results_all)

model_names = results_all.model_names;
contrast_labels = results_all.contrast_labels;

rows = {};

for m = 1:numel(model_names)
    model_name = model_names{m};
    res_model = results_all.(model_name);

    for c = 1:numel(contrast_labels)
        res = res_model{c};

        counts = accumarray(res.cluster_labels, 1);
        counts_str = mat2str(counts(:)');

        max_sil = max(res.mean_silhouette);

        rows(end+1, :) = { ...
            model_name, ...
            contrast_labels{c}, ...
            res.best_k, ...
            max_sil, ...
            counts_str ...
            };
    end
end

T = cell2table(rows, ...
    'VariableNames', {'Model', 'Contrast', 'BestK', 'MeanSilhouette', 'ClusterCounts'});

end