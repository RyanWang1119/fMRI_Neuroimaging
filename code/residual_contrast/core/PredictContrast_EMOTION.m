function PredictContrast_EMOTION(DATA_FILE, OUTDIR)
COMP2 = [1 2];                    % fear vs neut
cond_labels = {'fear','neut'};
run_observed_stage_flexible(DATA_FILE, OUTDIR, COMP2, cond_labels, 'EMOTION');
end
