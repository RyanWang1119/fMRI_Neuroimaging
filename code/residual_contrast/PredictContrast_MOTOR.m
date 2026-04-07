function PredictContrast_MOTOR(DATA_FILE, OUTDIR)
COMP2 = nchoosek(1:6,2);          % all pairwise (15)
cond_labels = {'cue','lf','rf','lh','rh','tongue'};
run_observed_stage_flexible(DATA_FILE, OUTDIR, COMP2, cond_labels, 'MOTOR');
end
