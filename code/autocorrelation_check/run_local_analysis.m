% This script uses a relative path to find the data, making it more robust.

% --- Make sure SPM is on your MATLAB path ---
% addpath('C:\path\to\your\spm12'); % <-- IMPORTANT: Update this path for Windows!

% --- Define Paths Relative to This Script ---
% Get the path to the folder containing this script ('.../code/')
scriptDir = fileparts(mfilename('fullpath'));

% Define the data directory as the 'data' folder next to the 'code' folder
dataDir = fullfile(scriptDir, '..', 'data');

% --- Add a check to make sure the data directory actually exists ---
if ~isfolder(dataDir)
    error('Data directory not found. Please ensure your project has this structure: Project/code and Project/data');
end

% Define a safe output location on the Desktop
outputBaseDir = fullfile(getenv('USERPROFILE'), 'Desktop', 'Autocorr_Outputs');

filesToProcess = {
    fullfile(dataDir, 'WM_100307_res4d.nii'),
    fullfile(dataDir, 'WM_108121_res4d.nii'),
    fullfile(dataDir, 'WM_200917_res4d.nii')
    };

% --- Loop through files and run the analysis ---
for i = 1:length(filesToProcess)
    analyze_local_subject_autocorr(filesToProcess{i}, outputBaseDir);
end

fprintf('All subjects processed! Check your Desktop for the Autocorr_Outputs folder.\n');