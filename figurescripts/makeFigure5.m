function makeFigure5(varargin)
% This is a function to make Figure 5 from the manuscript:
%   A visual encoding model links magnetoencephalography
%   signals to neural synchrony in human cortex.
%       by Kupers, Benson, Winawer (YEAR) JOURNAL.
%
% This figure shows the MEG forward model based on coherent and incoherent
% predictions coming from vertices located in V1-V3.
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
%   [plotMeanSubject]       :  (bool) true/false plot average across subjects
%   [saveFig]               :  (bool) true/false save figures
%   [headmodelType]         :  (str)  type of headmodel. Choose from 'OS' 
%                                     (overlapping spheres) or 'BEM'
%                                     (boundary element model)
%   [highResSurf]           :  (bool) true/false use high resolution 
%                                     headmodel/surface resolution.
%   [area]                  :  (str)  visual area to use from Benson et al.
%                                     (2014) PloS Comp Bio template. 
%                                     Choose from 'V1', 'V2', 'V3', 'V123'
%                                     'benson18atlas', or 'wang15atlas' (note:
%                                     eccentricity boundaries cannot be
%                                     defined when using wang15atlas)
%   [eccenLimitDeg]         :  (int)  eccentricity limit (deg) of the template. 
%                                     Supposingly matching the stimulus aperture. 
%                                     Can be a single int x, to get [0 x] or 
%                                     a vector [x,y] limiting eccentricity to
%                                     larger/equal to x and smaller/equal to y)
%   [contourPercentile]     :  (int)  percentile of the data to draw contours
%                                     If colormap max is at 97.5th percentile: 
%                                     Top15 sensors = 90.4. Top10 sensors = 93.6
%                                     To get contour lines at equal percentiles 
%                                     of data, use any integer under 10.
%   [maxColormapPercentile] :  (int)  percentile of data to truncate colormap
%   [signedColorbar]        :  (bool) true/false plot signed colormap or only
%                                     positive values.
%   [singleColorbarFlag]    :  (bool) Use a single colorbar per row,instead of 
%                                     per individual plot (default: true)
%
% Example 1:
%  makeFigure4('subjectsToPlot', 1, 'plotMeanSubject', false, 'saveFig', true)
% Example 2:
%  makeFigure4('subjectsToPlot', 12, 'plotMeanSubject', false, 'saveFig', true)
% Example 3:
%  makeFigure4('subjectsToPlot', 1:12, 'plotMeanSubject', true, 'saveFig', true)
%
% By Eline Kupers (NYU) 2017

%% 0. Set up paths and define parameters

p = inputParser;
p.KeepUnmatched = true;
p.addParameter('subjectsToPlot', 12);
p.addParameter('plotMeanSubject', true, @islogical)
p.addParameter('saveFig', true, @islogical);
p.addParameter('headmodelType', 'OS', @(x) any(validatestring(x,{'OS', 'BEM'})));
p.addParameter('highResSurf', false, @islogical);
p.addParameter('area', 'V123', @(x) any(validatestring(x,{'V1', 'V2', 'V3','V123', 'benson18atlas', 'wang15atlas'})));
p.addParameter('eccenLimitDeg', [0.18 11], @isnumeric);
p.addParameter('contourPercentile', 93.6, @isnumeric);
p.addParameter('maxColormapPercentile', 100, @isnumeric);
p.addParameter('signedColorbar', false, @islogical);
p.addParameter('singleColorbarFlag', true, @islogical);
p.parse(varargin{:});

% Rename variables
subjectsToPlot        = p.Results.subjectsToPlot;
plotMeanSubject       = p.Results.plotMeanSubject;
saveFig               = p.Results.saveFig;
headmodelType         = p.Results.headmodelType;
highResSurf           = p.Results.highResSurf;
area                  = p.Results.area;
eccenLimitDeg         = p.Results.eccenLimitDeg;
contourPercentile     = p.Results.contourPercentile;
maxColormapPercentile = p.Results.maxColormapPercentile;
signedColorbar        = p.Results.signedColorbar;
singleColorbarFlag    = p.Results.singleColorbarFlag;


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

% Load all subjects when plotting the mean
if plotMeanSubject
    subjectsToLoad = 1:length(subject);
else
    subjectsToLoad = subjectsToPlot;
end

% Plotting params
colorConds = {'y','b'};
markerType   = '.';


%% Loop over subjects
for s = subjectsToLoad
    
    d = dir(fullfile(bsDB, projectName, 'data', subject{s}, 'R*'));
    bsData = fullfile(d(1).folder, d(1).name);
    
    if highResSurf
        bsAnat = fullfile(bsDB, projectName, 'anat', subject{s}, 'highres');
    else
        bsAnat = fullfile(bsDB, projectName, 'anat', subject{s}, 'lowres');
    end
    %% 1. Load relevant matrices
    
    G_constrained = getGainMatrix(bsData, keep_sensors, headmodelType, highResSurf);
    
    % Get V1 template limited to 11 degrees eccentricity
    template = getTemplate(bsAnat, area, eccenLimitDeg);
    
    % Simulate coherent, in between or mixture, adn incoherent source time
    % series and compute predictions from forward model (w)
    tmp = getForwardModelPredictions(G_constrained, template.([area '_StimEccen']), [], n, nrEpochs, theta, kappa);
    
    % Compute amplitude across time
    amps.c = abs(fft(tmp.c,[],2));
    amps.i = abs(fft(tmp.i,[],2));
    
    % Compute mean weights across epochs at input frequency
    w.V123c(s,:) = mean(amps.c(:,2,:),3);
    w.V123i(s,:) = mean(amps.i(:,2,:),3);
    
end

%% Visualize predictions
    
% Take mean across subjects
w.V123c_mn = mean(w.V123c,1);
w.V123i_mn = mean(w.V123i,1);

% Get max of synchronous prediction to normalize all data to
maxSynAverage = max(w.V123c_mn);
    
% Plot average across subject if requested
if plotMeanSubject
    
    w.V123c_mn_norm = w.V123c_mn./maxSynAverage;
    w.V123i_mn_norm = w.V123i_mn./maxSynAverage;
    
    dataToPlot = cat(1,w.V123c_mn_norm, w.V123i_mn_norm);
    
    fig_ttl    = {sprintf('Figure5_ModelPredictions_%s_%1.2f-%d_prctile%d_%s_highResFlag%d_singleColorbarFlag%d_AVERAGE', ...
                        area, eccenLimitDeg(1),eccenLimitDeg(2), contourPercentile, headmodelType, highResSurf,singleColorbarFlag), ...
                  sprintf('Figure5_Contour_%s_%1.2f-%d_prctile%d_%s_highResFlag%d_singleColorbarFlag%d_AVERAGE', ...
                        area, eccenLimitDeg(1),eccenLimitDeg(2), contourPercentile, headmodelType, highResSurf,singleColorbarFlag)};
                    
    sub_ttl    = {sprintf('Synchronous sources - Average N = %d', length(subject)), ...
                  sprintf('Asynchronous sources - Average N = %d', length(subject))};
    
    figureDir  = fullfile(fmsRootPath,'figures', 'average'); % Where to save images?    
    if ~exist(figureDir,'dir'); mkdir(figureDir); end
    
    % Plot data and save channels that are located inside the contour lines
    visualizeSensormaps(dataToPlot, maxColormapPercentile, contourPercentile, signedColorbar, singleColorbarFlag, colorConds, markerType, fig_ttl, sub_ttl, saveFig, figureDir);
    
end

%% For individual subjects, normalized by mean
for s = subjectsToPlot
    synData = w.V123c(s,:)./maxSynAverage;
    asynData = w.V123i(s,:)./maxSynAverage;
    
    dataToPlot   = cat(1,synData, asynData);

    sub_ttl      = {sprintf('Synchronous sources S%d', s), ...
                    sprintf('Asynchronous sources S%d', s)};
                
    fig_ttl      = {sprintf('Figure5_ModelPredictions_%s_%1.2f-%d_prctile%2.1f_%s_highResFlag%d_singleColorbarFlag%d_S%d', ...
                        area, eccenLimitDeg(1),eccenLimitDeg(2), contourPercentile, headmodelType, highResSurf, singleColorbarFlag,s), ...
                    sprintf('Figure5_Contour_%s_%1.2f-%d_prctile%2.1f_%s_highResFlag%d_singleColorbarFlag%d_S%d', ...
                        area, eccenLimitDeg(1),eccenLimitDeg(2), contourPercentile, headmodelType, highResSurf, singleColorbarFlag,s)};
    
    % Make figure and data dir for subject, if non-existing
    figureDir    = fullfile(fmsRootPath,'figures', subject{s}); % Where to save images?
    if ~exist(figureDir,'dir'); mkdir(figureDir); end
    
    visualizeSensormaps(dataToPlot, maxColormapPercentile, contourPercentile, ...
                        signedColorbar, singleColorbarFlag, colorConds, markerType, fig_ttl, sub_ttl, saveFig, figureDir);
    
end

% Check toolbox versions
% out = ver;
% save(fullfile(fmsRootPath, 'toolboxVersions.mat'),'out')

