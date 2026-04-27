function PredictContrast_GAMBLING(DATA_FILE, OUTDIR)
COMP2 = [1 2; 1 3; 2 3];          % loss vs neut, loss vs win, neut vs win
cond_labels = {'loss','neut','win'};
run_observed_stage_flexible(DATA_FILE, OUTDIR, COMP2, cond_labels, 'GAMBLING');
end
