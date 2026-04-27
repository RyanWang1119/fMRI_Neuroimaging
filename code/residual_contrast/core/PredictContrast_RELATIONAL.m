function PredictContrast_RELATIONAL(DATA_FILE, OUTDIR)
COMP2 = nchoosek(1:3,2);          % error/match/relation → 3 pairwise
cond_labels = {'error','match','relation'};
run_observed_stage_flexible(DATA_FILE, OUTDIR, COMP2, cond_labels, 'RELATIONAL');
end
