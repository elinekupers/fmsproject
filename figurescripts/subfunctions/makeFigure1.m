function makeFigure1(varargin)
% This is a function to make the last panel of Figure 1 from the manuscript:
%   A visual encoding model links magnetoencephalography
%   signals to neural synchrony in human cortex.
%       by Kupers, Benson, Winawer (YEAR) JOURNAL.
%
% This figure shows the MEG forward model based on coherent
% simulated sine waves coming from vertices located in V1-V3.
%
% To runs this script, you need:
% (1) Access to the SSMEG folder in the brainstorm data base
%     (on the winawerlab server under '/Projects/MEG/brainstorm_db/'
% (2) MEG_utils and Fieldtrip toolbox added to the paths. For example:
%        tbUse('ForwardModelSynchrony');
%     or to only add the MEG_utils toolbox:
%        addpath(genpath('~/matlab/git/toolboxes/meg_utils'))
% (3) Run the s_visualAreasFS2BS script from this repository
%
% INPUTS:
%   [subjectsToPlot]        :  (int)  subject nr to plot, default is 12
%   [saveFig]               :  (bool) true/false save figures
%   [headmodelType]         :  (str)  type of headmodel. Choose from 'OS' 
%                                     (overlapping spheres) or 'BEM'
%                                     (boundary element model)
%   [highResSurf]           :  (bool) true/false use high resolution 
%                                     headmodel/surface resolution.
%   [area]                  :  (str)  visual area to use from Benson et al.
%                                     (2014) PloS Comp Bio template. 
%                                     Choose from 'V1', 'V2', 'V3', 'V123',
%                                     'V23', or
%                                     'benson18atlas', or 'wang15atlas' (note:
%                                     eccentricity boundaries cannot be
%                                     defined when using wang15atlas)
%   [eccenLimitDeg]         :  (int)  eccentricity limit (deg) of the template. 
%                                     Supposingly matching the stimulus aperture. 
%                                     Can be a single int x, to get [0 x] or 
%                                     a vector [x,y] limiting eccentricity to
%                                     larger/equal to x and smaller/equal to y)
%   [useConstrainedDipoles] :  (bool) Use constrained dipoles for gain
%                                     matrix/headmodel, i.e. perpendicular 
%                                     to surface (default = true) or
%                                     unconstrained diples (false) 
%
% Example 1:
%  makeFigure1('subjectsToPlot', 1, 'saveFig', true)
%
% By Eline Kupers (NYU) 2017

%% 0. Set up paths and define parameters

p = inputParser;
p.KeepUnmatched = true;
p.addParameter('subjectsToPlot', 2);
p.addParameter('saveFig', true, @islogical);
p.addParameter('headmodelType', 'OS', @(x) any(validatestring(x,{'OS', 'BEM'})));
p.addParameter('highResSurf', false, @islogical);
p.addParameter('area', 'V123', @(x) any(validatestring(x,{'V1', 'V2', 'V3','V23','V123', 'benson18atlas', 'wang15atlas'})));
p.addParameter('eccenLimitDeg', [0.18 11], @isnumeric);
p.addParameter('useConstrainedDipoles', true, @islogical);
p.parse(varargin{:});

% Rename variables
subjectsToPlot        = p.Results.subjectsToPlot;
saveFig               = p.Results.saveFig;
headmodelType         = p.Results.headmodelType;
highResSurf           = p.Results.highResSurf;
area                  = p.Results.area;
eccenLimitDeg         = p.Results.eccenLimitDeg;
useConstrainedDipoles = p.Results.useConstrainedDipoles;


%% 1. Define subjects and synchrony variables
% Get subject names and corresponding data session number
subject = getSubjectIDs;

% Path to brainstorm database and project name
bsDB            = '/Volumes/server/Projects/MEG/brainstorm_db/';
projectName     = 'SSMEG';

% Number of iterations for the random coherence prediction of the forward model
n        	= 10;        % number of timepoints (ms)
nrEpochs    = 1000;      % number of epochs
theta       = 0;         % von mises mean, equal for three distributions (syn, asyn and mix)
kappa.syn   = 100*pi;    % very narrow von Mises
kappa.asyn  = 0;         % very broad (uniform) von Mises
kappa.mix   = 0.27*pi;   % in-between width size von Mises (note: we are not plotting these values for this figure)

% Define vector that can truncate number of sensors
keep_sensors = logical([ones(157,1); zeros(192-157,1)]); % TODO: Figure out a more generic way to define keep_sensors
    
% Get brainstorm data
d = dir(fullfile(bsDB, projectName, 'data', subject{subjectsToPlot}, 'R*'));
bsData = fullfile(d(1).folder, d(1).name);
bsAnat = fullfile(bsDB, projectName, 'anat', subject{subjectsToPlot}, 'lowres');

%% 1. Load relevant matrices
G = getGainMatrix(bsData, keep_sensors, headmodelType, highResSurf, useConstrainedDipoles);

% Get V123 template limited to 11 degrees eccentricity
template = getTemplate(bsAnat, area, eccenLimitDeg);

% Simulate coherent, in between or mixture, adn incoherent source time
% series and compute predictions from forward model (w)
tmp = getForwardModelPredictions(G, template.([area '_StimEccen']), [], n, nrEpochs, theta, kappa);

% Get mean across epochs for several time points
for t = 1:n
    prediction(t,:) = squeeze(mean(tmp.c(:,t,:),3,'omitnan'));
end    


%% Visualize predictions
figure(1); clf; hold on;
for t = 1:n
    subplot(2,n/2,t);
    [~,ch] = megPlotMap(prediction(t,:),[-2 2].*10^-4,[],'bipolar',sprintf('t = %d',t),[],[], ...
        'mask', 'convex', ...
        'pointsize', 10);
    set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);
end


