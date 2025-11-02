tmap_fdr = fmri_data('data\group_t_map_fdr_05.nii');
orthviews(tmap_fdr);
regions_report = region(tmap_fdr);
orthviews(tmap_fdr); 
sgtitle('Population t map, FDR corrected');
saveas(gcf, 'Population_t-value_map_corrected.png');