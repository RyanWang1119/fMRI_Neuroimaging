function save_roi_influence_outputs(cfg, runDir, windowTRs, windowSeconds, subjectIDs, ...
    roiLabels, outerSplits, outerResults, summaryTable, signFlipResults, ...
    runMetadata, softwareVersionInfo, timestamp)
%SAVE_ROI_INFLUENCE_OUTPUTS Save MAT and CSV outputs for one model/contrast.

if ~exist(runDir, 'dir')
    mkdir(runDir);
end

contrastSafe = sanitize_name(runMetadata.contrastLabel);
matFile = fullfile(runDir, sprintf('roi_influence_%s_%s.mat', cfg.hrfModelName, contrastSafe));
csvFile = fullfile(runDir, 'summaryTable.csv');
readmeFile = fullfile(runDir, 'README_roi_influence.txt');

save(matFile, ...
    'cfg', 'windowTRs', 'windowSeconds', 'subjectIDs', 'roiLabels', ...
    'outerSplits', 'outerResults', 'summaryTable', 'signFlipResults', ...
    'runMetadata', 'softwareVersionInfo', 'timestamp', '-v7.3');

writetable(summaryTable, csvFile);
write_run_readme(readmeFile, cfg, runMetadata, matFile, csvFile);

fprintf('Saved MAT output: %s\n', matFile);
fprintf('Saved CSV summary: %s\n', csvFile);
end

function write_run_readme(readmeFile, cfg, meta, matFile, csvFile)
fid = fopen(readmeFile, 'w');
if fid < 0
    warning('Could not write run README: %s', readmeFile);
    return
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'WM ROI Influence Elastic-Net Output\n');
fprintf(fid, '===================================\n\n');
fprintf(fid, 'HRF model: %s\n', cfg.hrfModelName);
fprintf(fid, 'Contrast: %s\n', meta.contrastLabel);
fprintf(fid, 'Input file: %s\n', meta.dataFile);
fprintf(fid, 'MAT output: %s\n', matFile);
fprintf(fid, 'CSV summary: %s\n\n', csvFile);
fprintf(fid, 'Method summary:\n');
fprintf(fid, '- Subject-level paired differences are built as A=(2bk category A - 0bk category A), B=(2bk category B - 0bk category B), D=A-B.\n');
fprintf(fid, '- Elastic-net logistic regression uses [D/2; -D/2] with grouped subject-level inner folds and held-out subject-level TAFC.\n');
fprintf(fid, '- Haufe pattern, elastic-net selection frequency, and held-out model reliance are separate outputs.\n');
fprintf(fid, '- The sign-flip maxT test is independent of decoder coefficients and selections.\n\n');
fprintf(fid, 'Interpretation:\n');
fprintf(fid, 'The high_confidence_influential_candidate flag is a reproducibility-based candidate label.\n');
fprintf(fid, 'It is not a causal claim and should not be read as mechanistic, exclusive, or definitive biological evidence.\n');
fprintf(fid, 'Held-out model reliance is not a causal measure; correlated parcels can substitute for each other.\n');
end

function s = sanitize_name(s)
s = char(s);
s = regexprep(s, '[^A-Za-z0-9_+-]+', '_');
s = regexprep(s, '_+', '_');
s = regexprep(s, '^_|_$', '');
end
