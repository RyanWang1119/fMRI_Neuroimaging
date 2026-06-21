function outerSplits = make_repeated_subject_splits(nSubjects, subjectIDs, cfg, splitFile)
%MAKE_REPEATED_SUBJECT_SPLITS Repeated subject-level subsampling splits.
%
% Splits are saved once and reused across HRF models and interaction
% contrasts. The unit of splitting is always the original subject.

if nargin < 4 || isempty(splitFile)
    splitFile = '';
end

subjectIDs = subjectIDs(:);
if numel(subjectIDs) ~= nSubjects
    error('subjectIDs length (%d) must match nSubjects (%d).', numel(subjectIDs), nSubjects);
end

if ~isempty(splitFile) && isfile(splitFile)
    S = load(splitFile, 'outerSplits', 'splitMetadata');
    if ~isfield(S, 'outerSplits') || ~isfield(S, 'splitMetadata')
        error('Shared split file is missing outerSplits or splitMetadata: %s', splitFile);
    end
    validate_split_metadata(S.splitMetadata, nSubjects, subjectIDs, cfg);
    outerSplits = S.outerSplits;
    validate_outer_splits(outerSplits, nSubjects);
    return
end

nRepeats = cfg.outer.nRepeats;
trainFraction = cfg.outer.trainFraction;
if trainFraction <= 0 || trainFraction >= 1
    error('cfg.outer.trainFraction must be between 0 and 1.');
end

nTrain = floor(trainFraction * nSubjects);
if nTrain < 2 || nSubjects - nTrain < 1
    error('Invalid train/test sizes for nSubjects=%d trainFraction=%g.', nSubjects, trainFraction);
end

outerSplits = repmat(struct( ...
    'repeat', NaN, ...
    'seed', NaN, ...
    'trainIdx', [], ...
    'testIdx', [], ...
    'trainSubjectIDs', [], ...
    'testSubjectIDs', []), nRepeats, 1);

for r = 1:nRepeats
    repeatSeed = cfg.seed + r * 1009;
    stream = RandStream('mt19937ar', 'Seed', repeatSeed);
    perm = randperm(stream, nSubjects);
    trainIdx = sort(perm(1:nTrain))';
    testIdx = sort(perm(nTrain+1:end))';

    outerSplits(r).repeat = r;
    outerSplits(r).seed = repeatSeed;
    outerSplits(r).trainIdx = trainIdx;
    outerSplits(r).testIdx = testIdx;
    outerSplits(r).trainSubjectIDs = subjectIDs(trainIdx);
    outerSplits(r).testSubjectIDs = subjectIDs(testIdx);
end

validate_outer_splits(outerSplits, nSubjects);

splitMetadata = struct();
splitMetadata.created = datestr(now, 30);
splitMetadata.seed = cfg.seed;
splitMetadata.nRepeats = nRepeats;
splitMetadata.trainFraction = trainFraction;
splitMetadata.nSubjects = nSubjects;
splitMetadata.subjectIDs = subjectIDs;
splitMetadata.note = 'Subject-level splits reused across HRF-model and interaction-contrast runs.';

if ~isempty(splitFile)
    splitDir = fileparts(splitFile);
    if ~exist(splitDir, 'dir')
        mkdir(splitDir);
    end
    save(splitFile, 'outerSplits', 'splitMetadata', '-v7.3');
end
end

function validate_split_metadata(meta, nSubjects, subjectIDs, cfg)
if meta.nSubjects ~= nSubjects
    error('Shared split file has nSubjects=%d, current data has nSubjects=%d.', meta.nSubjects, nSubjects);
end
if meta.nRepeats ~= cfg.outer.nRepeats
    error('Shared split file has nRepeats=%d, cfg requests %d.', meta.nRepeats, cfg.outer.nRepeats);
end
if abs(meta.trainFraction - cfg.outer.trainFraction) > eps
    error('Shared split file trainFraction=%g, cfg requests %g.', meta.trainFraction, cfg.outer.trainFraction);
end
if meta.seed ~= cfg.seed
    error('Shared split file seed=%d, cfg requests %d.', meta.seed, cfg.seed);
end
if numel(meta.subjectIDs) ~= numel(subjectIDs) || ~isequal(meta.subjectIDs(:), subjectIDs(:))
    error('Shared split file subjectIDs do not match current data. Use the same subject order or remove the split file intentionally.');
end
end

function validate_outer_splits(outerSplits, nSubjects)
for r = 1:numel(outerSplits)
    tr = outerSplits(r).trainIdx(:);
    te = outerSplits(r).testIdx(:);
    if any(tr < 1 | tr > nSubjects) || any(te < 1 | te > nSubjects)
        error('Outer split %d contains subject indices outside 1..%d.', r, nSubjects);
    end
    if numel(unique(tr)) ~= numel(tr) || numel(unique(te)) ~= numel(te)
        error('Outer split %d contains duplicate subject indices.', r);
    end
    if ~isempty(intersect(tr, te))
        error('Outer split %d has subjects in both train and test sets.', r);
    end
end
end
