function makeFigure8(varargin)
% This is a function to make Figure 8 from the manuscript:
%   A visual encoding model links magnetoencephalography
%   signals to neural synchrony in human cortex.
%       by Kupers, Benson, Winawer (YEAR) JOURNAL.
%
% This figure shows the MEG forward model based on coherent and incoherent
% predictions coming from vertices located in V1 when using an ABSOLUTE
% gain matrix (so propagation of magnetic fields WITHOUT CANCELLATION).
%
% To runs this script, you need: (1) Access to the SSMEG folder in the
% brainstorm data base (2) MEG_utils and Fieldtrip toolbox added to the
% paths. For example:
%     tbUse('ForwardModelSynchrony');
%        or to only add the MEG_utils toolbox:
%     addpath(genpath('~/matlab/git/toolboxes/meg_utils'))
%
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
%                                     'benson18', or 'wang15atlas' (note:
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
%   [singleColorbarFlag]    :  (bool) Use a single colorbar per row,
%                                     instead of per individual plot.
% Example 1:
%  makeFigure8('subjectsToPlot', 1, 'plotMeanSubject', false, 'saveFig', true)
% Example 2:
%  makeFigure8('subjectsToPlot', 12, 'plotMeanSubject', false, 'saveFig', true)
% Example 3:
%  makeFigure8('subjectsToPlot', 1:12, 'plotMeanSubject', true, 'saveFig', true)
%
%
% By Eline Kupers (NYU) 2019

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
p.addParameter('useConstrainedDipoles', true, @islogical);
p.addParameter('contourPercentile', 93.6, @isnumeric);
p.addParameter('maxColormapPercentile', 100, @isnumeric);
p.addParameter('signedColorbar', false, @islogical);
p.addParameter('singleColorbarFlag', false, @islogical);
p.parse(varargin{:});

% Rename variables
subjectsToPlot        = p.Results.subjectsToPlot;
plotMeanSubject       = p.Results.plotMeanSubject;
saveFig               = p.Results.saveFig;
headmodelType         = p.Results.headmodelType;
highResSurf           = p.Results.highResSurf;
area                  = p.Results.area;
eccenLimitDeg         = p.Results.eccenLimitDeg;
useConstrainedDipoles = p.Results.useConstrainedDipoles;
contourPercentile     = p.Results.contourPercentile;
maxColormapPercentile = p.Results.maxColormapPercentile;
signedColorbar        = p.Results.signedColorbar;
singleColorbarFlag    = p.Results.singleColorbarFlag;


% Get subject names
subject = getSubjectIDs;

% Load all subjects when plotting the mean
if plotMeanSubject
    subjectsToLoad = 1:length(subject);
else
    subjectsToLoad = subjectsToPlot;
end

% Path to brainstorm database and project name
bsDB            = '/Volumes/server/Projects/MEG/brainstorm_db/';
projectName     = 'SSMEG';

% Define params for forward model
n        	= 10;        % number of timepoints (ms)
nrEpochs    = 1000;      % number of epochs
theta       = 0;         % von mises mean, equal for three distributions
kappa.syn   = 100*pi;    % very narrow von Mises
kappa.asyn  = 0;         % very broad (uniform) von Mises
kappa.mix   = 0.27*pi;   % in between width size von Mises

% Define vector that can truncate number of sensors
keep_sensors = logical([ones(157,1); zeros(192-157,1)]); % TODO: Figure out a more generic way to define keep_sensors

% Define plotting params
colorMarkers = {'y','b'};
markerType   = '.';

%% Loop over subjects
for s = subjectsToLoad
    
    % Get subject data and anatomy files
    d = dir(fullfile(bsDB, projectName, 'data', subject{s}, 'R*'));
    bsData = fullfile(d(1).folder, d(1).name);
    
    if highResSurf
        bsAnat = fullfile(bsDB, projectName, 'anat', subject{s}, 'highres');
    else
        bsAnat = fullfile(bsDB, projectName, 'anat', subject{s}, 'lowres');
    end
    
    %% 1. Load relevant matrices
    
    % Simulate coherent and incoherent source time series and compute
    % predictions from forward model (w)
    G_constrained = getGainMatrix(bsData, keep_sensors, headmodelType, highResSurf, useConstrainedDipoles);
    
    % Get V1 template limited to 11 degrees eccentricity
    template = getTemplate(bsAnat, area, eccenLimitDeg);
    
    % Simulate coherent, in between or mixture, adn incoherent source time
    % series and compute predictions from forward model (w)
    % IMPORTANT: We when we take absolute values of G_contrained there is no cancellation possible
    tmp.woC = getForwardModelPredictions(abs(G_constrained), template.([area '_StimEccen']), [], n, nrEpochs, theta, kappa);
    tmp.wC = getForwardModelPredictions(G_constrained, template.([area '_StimEccen']), [], n, nrEpochs, theta, kappa);
    
    % Compute amplitude at freq
    amps.woC.c = abs(fft(tmp.woC.c,[],2));
    amps.woC.i = abs(fft(tmp.woC.i,[],2));
    
    amps.wC.c = abs(fft(tmp.wC.c,[],2));
    amps.wC.i = abs(fft(tmp.wC.i,[],2));
    
    % Take mean across epochs
    w.woC.V123c(s,:) = mean(amps.woC.c(:,2,:),3);
    w.woC.V123i(s,:) = mean(amps.woC.i(:,2,:),3);
    
    w.wC.V123c(s,:) = mean(amps.wC.c(:,2,:),3);
    w.wC.V123i(s,:) = mean(amps.wC.i(:,2,:),3);
    
end

%% 2. Get average across subjects

% Take the average across subjects (with Cancellation)
w.wC.V123c_mn       = mean(w.wC.V123c,1);
w.wC.V123i_mn       = mean(w.wC.V123i,1);

maxSynAverage = max(w.wC.V123c_mn);
w.wC.V123c_mn_norm  = w.wC.V123c_mn./maxSynAverage;
w.wC.V123i_mn_norm  = w.wC.V123i_mn./maxSynAverage;

% Take the average across subjects (without Cancellation)
w.woC.V123c_mn      = mean(w.woC.V123c,1);
w.woC.V123i_mn      = mean(w.woC.V123i,1);

w.woC.V123c_mn_norm = w.woC.V123c_mn./maxSynAverage;
w.woC.V123i_mn_norm = w.woC.V123i_mn./maxSynAverage;


% Take ratio
ratioCoh_norm   = w.woC.V123c_mn_norm ./ w.wC.V123c_mn_norm;
ratioInCoh_norm = w.woC.V123i_mn_norm ./ w.wC.V123i_mn_norm;

% Visualize predictions from forward model for group average
if plotMeanSubject
    
    % Redefine figure dir
    figureDir    = fullfile(fmsRootPath, 'figures', 'average'); % Where to save images?
    
    % Define plotting data and figure titles
    dataToPlot   = cat(1, w.woC.V123c_mn_norm, w.woC.V123i_mn_norm);
    
    fig_ttl      = {sprintf('Figure8_ModelPredictions-No_cancellation_%s_%1.2f-%d_prctle%2.1f_%s_highResFlag%d_singleColorbarFlag%d_AVERAGE', area, eccenLimitDeg(1),eccenLimitDeg(2), contourPercentile, headmodelType, highResSurf,singleColorbarFlag), ...
                    sprintf('Figure8_Contours-No_cancellation_%s_%1.2f-%d_prctle%2.1f_%s_highResFlag%d_singleColorbarFlag%d_AVERAGE', area, eccenLimitDeg(1),eccenLimitDeg(2), contourPercentile, headmodelType, highResSurf,singleColorbarFlag)};
    sub_ttl      = {sprintf('No cancellation: Synchronous sources Average N=%d', length(subject)), ...
                    sprintf('No cancellation: Asynchronous sources Average N=%d', length(subject))};
    
    % Plot it!
    visualizeSensormaps(dataToPlot, maxColormapPercentile, contourPercentile, ...
        signedColorbar, singleColorbarFlag, colorMarkers, markerType, fig_ttl, sub_ttl, saveFig, figureDir);
    
    % Plot log of ratio of with vs without cancellation
    dataToPlot = cat(1, ratioCoh_norm, ratioInCoh_norm);
    fig_ttl      = {sprintf('Figure8_ratioWithVsWithoutCancellation_%s_%1.2f-%d_prctile%2.1f_%s_highResFlag%d_singleColorbarFlag%d_AVERAGE', area, eccenLimitDeg(1),eccenLimitDeg(2), contourPercentile, headmodelType, highResSurf,singleColorbarFlag)};
    sub_ttl      = {'Syn: ratio with / without', ...
                    'Asyn: ratio with / without'};
  
    visualizeSensormaps(dataToPlot, 100, [], signedColorbar, singleColorbarFlag, colorMarkers, ...
            markerType, fig_ttl, sub_ttl, saveFig, figureDir);
        
    % Plot 1D average of posterior sensors 
    for s = 1:size(w.woC.V123c,1)
        synDataWoC = w.woC.V123c(s,:)./maxSynAverage;
        asynDataWoC = w.woC.V123i(s,:)./maxSynAverage;
        synDataWC = w.wC.V123c(s,:)./maxSynAverage;
        asynDataWC = w.wC.V123i(s,:)./maxSynAverage;
        
        allDataWoC{s} = cat(1,synDataWoC, asynDataWoC);
        allDataWC{s} = cat(1,synDataWC, asynDataWC);

    end
        
    fig_ttl2 = sprintf('Figure8_1Daverage_ModelPredictionsWithoutCancellation_AVERAGE');
    sub_ttl      = {'Syn: without cancellation', ...
                    'Asyn: without cancellation'};
    visualizePosteriorSensors1D(allDataWoC, plotMeanSubject, fig_ttl2, sub_ttl, saveFig, figureDir)

    fig_ttl3 = sprintf('Figure8_1Daverage_ModelPredictionsWithCancellation_AVERAGE');
    sub_ttl      = {'Syn: with cancellation', ...
                    'Asyn: with cancellation'};
    visualizePosteriorSensors1D(allDataWC, plotMeanSubject, fig_ttl3, sub_ttl, saveFig, figureDir) 
    
end

%% 3. Plot for requested individual subject
for exampleSubject = subjectsToPlot
    
    dataToPlot   = cat(1, (w.woC.V123c(exampleSubject,:)./maxSynAverage), ...
                          (w.woC.V123i(exampleSubject,:)./maxSynAverage));
    
    fig_ttl      = {sprintf('Figure8_ModelPredictions-No_cancellation_%s_%1.2f-%d_prctle%d_%s_highResFlag%d_singleColorbarFlag%d_S%d', ...
                      area, eccenLimitDeg(1),eccenLimitDeg(2), contourPercentile,headmodelType, highResSurf, singleColorbarFlag, exampleSubject), ...
                    sprintf('Figure8_Contours-No_cancellation_%s_%1.2f-%d_prctle%d_%s_highResFlag%d_singleColorbarFlag%d_S%d', ...
                      area, eccenLimitDeg(1),eccenLimitDeg(2), contourPercentile, headmodelType, highResSurf, singleColorbarFlag,exampleSubject)};
    sub_ttl      = {sprintf('No cancellation: Synchronous sources S%d', exampleSubject), ...
                    sprintf('No cancellation: Asynchronous sources S%d', exampleSubject)};
        
    % Make figure and data dir for subject, if non-existing
    figureDir    = fullfile(fmsRootPath,'figures', subject{exampleSubject}); % Where to save images?
    if ~exist(figureDir,'dir'); mkdir(figureDir); end
    
    % Plot it!
    visualizeSensormaps(dataToPlot, maxColormapPercentile, contourPercentile, ...
        signedColorbar, singleColorbarFlag, colorMarkers, markerType, fig_ttl, sub_ttl, saveFig, figureDir);
end


return