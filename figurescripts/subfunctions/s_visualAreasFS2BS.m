%% s_visualAreasFS2BS

% 0. Before you run this script:
% 0.1 Apply Freesurfer's recon-all all to subject's T1 anatomy
% 0.2 Run Benson template V123 with Docker on Freesurfer subject
% 0.3 Add subject to Brainstorm database (Freesurfer folder to anatomy) and
%    MEG data to functional data tab.
% 0.4 Add fiducials to T1
% 0.5 Check MRI/MEG registration

% 1. Flow of the script:
%
% 1.0 Downsample Benson atlas, to Brainstorm surface size. And use downsampled
% template to create a surface only containing V1,2,3 (unitized).
% 1.1 Export the unitized surface to a mgz file


% Toolbox dependencies:
% - Add Brainstorm, meg_utils, fieltrip toolbox by typing tbUse('ForwardModelSynchrony');
% and freesurfer matlab functions by :
%   addpath(genpath('/Applications/freesurfer/matlab'));
%   addpath(genpath('/Applications/freesurfer/fsfast/toolbox'))
% - Have the brainstorm GUI open


%% 0. Define paths and variables
% subject         = {'wlsubj002', ... % S1 - Full, Left, Right stim experiment
%     'wlsubj004', ... % S2 - Full, Left, Right stim experiment
%     'wlsubj005', ... % S3 - Full, Left, Right stim experiment
%     'wlsubj006', ... % S4 - Full, Left, Right stim experiment
%     'wlsubj010', ... % S5 - Full, Left, Right stim experiment
%     'wlsubj011', ... % S6 - Full, Left, Right stim experiment
%     'wlsubj048', ... % S7 - Full  stim only experiment
%     'wlsubj046', ... % S8 - Full  stim only experiment
%     'wlsubj039', ... % S9 - Full  stim only experiment
%     'wlsubj059', ... % S10 - Full  stim only experiment
%     'wlsubj067', ... % S11 - Full  stim only experiment
%     'wlsubj070'};    % S12 - Full  stim only experiment
% 

if ~exist('MRIread')
    addpath(genpath('/Applications/freesurfer/matlab'));
    addpath(genpath('/Applications/freesurfer/fsfast/toolbox'))
end

bssubject = 'wlsubj002';
fssubject = bssubject;
exp       = 'SSMEG';
highResFlag = false;

% Freesurfer, brainstorm database and data/anatomy directories
fsdir     = '/Volumes/server/Freesurfer_subjects';
bsDB     = '/Volumes/server/Projects/MEG/brainstorm_db';

d = dir(fullfile(bsDB, exp, 'data', bssubject));
datadir = fullfile(bsDB, exp, 'data', bssubject, d(end).name);

anatdir = fullfile(bsDB, exp, 'anat', bssubject);

%% 1.1 Create downsampled amplitude template
interp_retinotopy(bsDB, fsdir, fssubject, bssubject, exp, highResFlag)

