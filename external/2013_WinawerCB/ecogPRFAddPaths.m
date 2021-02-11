
function ecogPRFAddPaths
% Set Matlab directory path for ECOG PRF project
%
%     ecogPRFAddPaths
%
% Set up the path to the functions called by ECOG PRF functions.

rootPath = ecogPRFrootPath;
% fprintf('ecogPRFpath root directory: %s\n',rootPath)

% Adds the root directory of the ECOG PRF tree to the user's path
addpath(rootPath);

% Generates a list of the directories below the ECOG PRF tree.
addpath(fullfile(rootPath, 'data'))
addpath(genpath(fullfile(rootPath, 'external')))
addpath(fullfile(rootPath, 'figures'))
addpath(fullfile(rootPath, 'scratch'))
addpath(fullfile(rootPath, 'compute'))
addpath(fullfile(rootPath, 'simulation'))


return;

