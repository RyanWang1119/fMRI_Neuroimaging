function Tminor = get_minority_regions_from_result(res, atlas_obj)

counts = accumarray(res.cluster_labels, 1);
[~, minority_cluster] = min(counts);

minor_local_idx = find(res.cluster_labels == minority_cluster);
minor_orig_idx  = res.keep_idx(minor_local_idx);

if nargin > 1 && ~isempty(atlas_obj)
    minor_names = atlas_obj.labels(minor_orig_idx);
else
    if isfield(res, 'region_names') && ~isempty(res.region_names)
        minor_names = res.region_names(minor_local_idx);
    else
        minor_names = repmat({''}, numel(minor_orig_idx), 1);
    end
end

Tminor = table( ...
    repmat(minority_cluster, numel(minor_orig_idx), 1), ...
    minor_orig_idx(:), ...
    minor_names(:), ...
    'VariableNames', {'Cluster', 'ParcelIndex', 'ParcelName'});

end