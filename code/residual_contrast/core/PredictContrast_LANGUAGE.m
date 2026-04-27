function PredictContrast_LANGUAGE(DATA_FILE, OUTDIR)

COMP2 = nchoosek(1:6,2);          % all pairwise (15)
cond_labels = {'present_math','question_math','response_math','present_story','question_story','response_story'};
run_observed_stage_flexible(DATA_FILE, OUTDIR, COMP2, cond_labels, 'LANGUAGE');
end
