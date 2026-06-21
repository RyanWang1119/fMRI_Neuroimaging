function foldID = make_grouped_inner_folds(nTrainSubjects, innerKFolds, seed)
%MAKE_GROUPED_INNER_FOLDS Assign each training subject to one inner fold.
%
% The returned foldID is [nTrainSubjects x 1]. Downstream paired rows are
% made as [D/2; -D/2], so row membership is always derived from this subject
% fold vector and both synthetic rows for a subject remain together.

if nTrainSubjects < 2
    error('Need at least two training subjects for grouped inner folds.');
end
if innerKFolds < 2
    error('innerKFolds must be at least 2.');
end
innerKFolds = min(innerKFolds, nTrainSubjects);

stream = RandStream('mt19937ar', 'Seed', seed);
perm = randperm(stream, nTrainSubjects);
foldID = zeros(nTrainSubjects, 1);
for i = 1:nTrainSubjects
    foldID(perm(i)) = mod(i - 1, innerKFolds) + 1;
end

for k = 1:innerKFolds
    if ~any(foldID == k)
        error('Inner fold %d is empty.', k);
    end
end
end
