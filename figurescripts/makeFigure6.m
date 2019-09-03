function makeFigure6(varargin)
% This is a function to make Figure 6 from the manuscript about forward
% modeling coherent and incoherent neural sources to MEG responses.
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
%   [subjectsToPlot]  :  (int) subject nr you would like to plot, default is 12
%   [plotMeanSubject] :  (bool) plot average across all 12 subjets or not?
%   [saveFig]         :  (bool) save figures or not?
%
% Example 1:
%  makeFigure6('subjectsToPlot', 1, 'plotMeanSubject', false, 'saveFig', true)
% Example 2:
%  makeFigure6('subjectsToPlot', 12, 'plotMeanSubject', false, 'saveFig', true)
% Example 3:
%  makeFigure6('subjectsToPlot', 1:12, 'plotMeanSubject', true, 'saveFig', true)

%% 0. Set up paths and define parameters

p = inputParser;
p.KeepUnmatched = true;
p.addParameter('subjectsToPlot', 12);
p.addParameter('plotMeanSubject', true, @islogical)
p.addParameter('saveFig', true, @islogical);
p.addParameter('highResSurf', false, @islogical);                       % What surface resolution? Choose from true (high) or false (low)
p.addParameter('area', 'V123', @(x) any(x,{'V1', 'V2', 'V3','V123'}));  % What visual area to use? Choose between 'V1', 'V2', 'V3', or 'V123'
p.addParameter('eccenLimitDeg', [0.18 11], @isnumeric);                 % What is the eccentricity limit (deg) for the template? (supposingly matching the stimulus aperture.) 
                                                                        %   Can be a single int x, to get [0 x] or a vector [x,y] limiting eccentricity to larger/equal to x and smaller/equal to y)
p.addParameter('contourPercentile', 93.6, @isnumeric);                  % At what percentile of the data to draw contour lines?
                                                                        %   If colormap max is at 97.5th percentile: Top 15 channels -> 90.4. Top 10 channels -> 93.6,
                                                                        %   to get contour lines at equal percentiles of data, use any integer under 10
p.addParameter('maxColormapPercentile', 97.5, @isnumeric);              % At what percentile of data are we truncating colormap?
p.addParameter('signedColorbar', false, @islogical);                    % Plot signed colormap (true) or only positive values (false)? 
p.parse(varargin{:});

% Rename variables
subjectsToPlot        = p.Results.subjectsToPlot;
plotMeanSubject       = p.Results.plotMeanSubject;
saveFig               = p.Results.saveFig;
highResSurf           = p.Results.highResSurf;
area                  = p.Results.area;
eccenLimitDeg         = p.Results.eccenLimitDeg;
contourPercentile     = p.Results.contourPercentile;
maxColormapPercentile = p.Results.maxColormapPercentile;
signedColorbar        = p.Results.signedColorbar;


%% 1. Define subjects and synchrony variables
subject         = {'wlsubj002', ... S1 - From exp: Full, Left, Right
    'wlsubj004', ... S2 - From exp: Full, Left, Right
    'wlsubj005', ... S3 - From exp: Full, Left, Right
    'wlsubj006', ... S4 - From exp: Full, Left, Right
    'wlsubj010', ... S5 - From exp: Full, Left, Right
    'wlsubj011', ... S6 - From exp: Full, Left, Right
    'wlsubj048', ... S7 - From exp: Full only
    'wlsubj046', ... S8 - From exp: Full only
    'wlsubj039', ... S9 - From exp: Full only
    'wlsubj059', ... S10 - From exp: Full only
    'wlsubj067', ... S11 - From exp: Full only
    'wlsubj070'}; %  S12 - From exp: Full only

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

%% Loop over subjects
for s = subjectsToLoad
    
    d = dir(fullfile(bsDB, projectName, 'data', subject{s}, 'R*'));
    bsData = fullfile(d(1).folder, d(1).name);
    
    if highResSurf
        bsAnat = fullfile(bsDB, projectName, 'anat', subject{s}, 'highres');
    else
        bsAnat = fullfile(bsDB, projectName, 'anat', subject{s});
    end
    %% 1. Load relevant matrices
    
    G_constrained = getGainMatrix(bsData, keep_sensors, highResSurf);
    
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
colorConds = {'y','b'};
markerType   = '.';

for s = subjectsToPlot
    dataToPlot   = cat(1,w.V123c(s,:), w.V123i(s,:));

    sub_ttl      = {sprintf('Synchronous sources S%d', s), ...
                    sprintf('Asynchronous sources S%d', s)};
                
    fig_ttl      = {sprintf('Figure6_ModelPredictions_%s_%1.2f-%d_prctile%2.1f_highResFlag%d_S%d', ...
                        area, eccenLimitDeg(1),eccenLimitDeg(2), contourPercentile, highResSurf, s), ...
                    sprintf('Figure6_Contour_%s_%1.2f-%d_prctile%2.1f_highResFlag%d_S%d', ...
                        area, eccenLimitDeg(1),eccenLimitDeg(2), contourPercentile, highResSurf, s)};

    dataDir      = fullfile(fmsRootPath,'data', subject{s}); % Where to save vector of sensors that fall within contours?
    figureDir    = fullfile(fmsRootPath,'figures', subject{s}); % Where to save images?
    
    % Make figure and data dir for subject, if non-existing
    if ~exist(figureDir,'dir'); mkdir(figureDir); end
    if ~exist(dataDir,'dir'); mkdir(dataDir); end
    
    visualizeSensormaps(dataToPlot, maxColormapPercentile, contourPercentile, ...
                        signedColorbar, colorConds, markerType, fig_ttl, sub_ttl, saveFig, figureDir);
   
    
    % Save sensors of interest falling within the contour lines
    if saveFig
        save(fullfile(dataDir, ...
            sprintf('%s_sensorsWithinContours_Prediction_%s_%1.2f-%d_prctile%2.1f_highResFlag%d_S%d.mat', ...
            subject{s}, area, eccenLimitDeg(1),eccenLimitDeg(2), contourPercentile, highResSurf, s)), 'dataToPlot'); 
    end
    
end
%% Take mean across subjects and plot if requested
if plotMeanSubject
    
    w.V123c_mn = mean(w.V123c,1);
    w.V123i_mn = mean(w.V123i,1);
    
    dataToPlot = cat(1,w.V123c_mn, w.V123i_mn);
    
    fig_ttl    = {sprintf('Figure6_ModelPredictions_%s_%1.2f-%d_prctile%d_highResFlag%d_AVERAGE', ...
                        area, eccenLimitDeg(1),eccenLimitDeg(2), contourPercentile, highResSurf), ...
                  sprintf('Figure6_Contour_%s_%1.2f-%d_prctile%d_highResFlag%d_AVERAGE', ...
                        area, eccenLimitDeg(1),eccenLimitDeg(2), contourPercentile, highResSurf)};
                    
    sub_ttl    = {sprintf('Synchronous sources - Average N = %d', length(subject)), ...
                  sprintf('Asynchronous sources - Average N = %d', length(subject))};
    
    figureDir  = fullfile(fmsRootPath,'figures', 'average'); % Where to save images?
    dataDir    = fullfile(fmsRootPath,'data', 'average'); % Where to save images?
    
    if ~exist(figureDir,'dir'); mkdir(figureDir); end
    if ~exist(dataDir,'dir'); mkdir(dataDir); end
    
    % Plot data and save channels that are located inside the contour lines
    visualizeSensormaps(dataToPlot, maxColormapPercentile, contourPercentile, signedColorbar, colorConds, markerType, fig_ttl, sub_ttl, saveFig, figureDir);
    
    % Save sensors of interest falling within the contour lines
    if saveFig
        save(fullfile(dataDir, sprintf('sensorsWithinContours_Prediction_%s_%1.2f-%d_prctile%2.1f_highResFlag%d_AVERAGE.mat', ...
            area,eccenLimitDeg(1),eccenLimitDeg(2),contourPercentile,highResSurf)), 'dataToPlot');
    end
    
end

% Check toolbox versions
% out = ver;
% save(fullfile(fmsRootPath, 'toolboxVersions.mat'),'out')

