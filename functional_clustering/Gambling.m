atlas_obj = load_atlas('canlab2018');

%% GAMBLING_loss
GAMBLING_loss = permute(Results(:,:,:,1), [3 1 2]);   
loss = region_cluster(GAMBLING_loss);
subj_loss = subject_cluster(GAMBLING_loss);

counts_loss = accumarray(loss.cluster_labels, 1);
[~, minority_cluster_loss] = min(counts_loss);
%% loss
minor_local_idx_loss = find(loss.cluster_labels == minority_cluster_loss);
minor_orig_idx_loss  = loss.keep_idx(minor_local_idx_loss);

minor_names_loss = atlas_obj.labels(minor_orig_idx_loss);

table(minor_orig_idx_loss(:), minor_names_loss(:), ...
    'VariableNames', {'ParcelIndex', 'ParcelName'})

loss_subset = select_atlas_subset(atlas_obj, minor_orig_idx_loss);
loss_region = atlas2region(loss_subset);

figure;
orthviews(loss_region);

%% GAMBLING_win
GAMBLING_win = permute(Results(:,:,:,3), [3 1 2]);   
win = region_cluster(GAMBLING_win);
subj_win = subject_cluster(GAMBLING_win);

%% win
counts_win = accumarray(win.cluster_labels, 1);
[~, minority_cluster_win] = min(counts_win);

minor_local_idx_win = find(win.cluster_labels == minority_cluster_win);
minor_orig_idx_win  = win.keep_idx(minor_local_idx_win);

minor_names_win = atlas_obj.labels(minor_orig_idx_win);

table(minor_orig_idx_win(:), minor_names_win(:), ...
    'VariableNames', {'ParcelIndex', 'ParcelName'})

win_subset = select_atlas_subset(atlas_obj, minor_orig_idx_win);
win_region = atlas2region(win_subset);

figure;
orthviews(win_region);

%% GAMBLING_neutral
GAMBLING_neutral = permute(Results(:,:,:,2), [3 1 2]);   
neutral = region_cluster(GAMBLING_neutral);
subj_neutral = subject_cluster(GAMBLING_neutral);

%% neutral
counts_neutral = accumarray(neutral.cluster_labels, 1);
[~, minority_cluster_neutral] = min(counts_neutral);

minor_local_idx_neutral = find(neutral.cluster_labels == minority_cluster_neutral);
minor_orig_idx_neutral  = neutral.keep_idx(minor_local_idx_neutral);

minor_names_neutral = atlas_obj.labels(minor_orig_idx_neutral);

table(minor_orig_idx_neutral(:), minor_names_neutral(:), ...
    'VariableNames', {'ParcelIndex', 'ParcelName'})

neutral_subset = select_atlas_subset(atlas_obj, minor_orig_idx_neutral);
neutral_region = atlas2region(neutral_subset);

figure;
orthviews(neutral_region);