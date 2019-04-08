function makeFigure1(varargin)
% This is a function to make Figure 1 from the manuscript about forward
% modeling coherent and incoherent neural sources to MEG responses.
%
% This figure shows the MEG forward model based on coherent and incoherent
% predictions coming from vertices located in V1.
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
%   [exampleSubject]  :  (int) subject nr you would like to plot, default is 12
%   [plotMeanSubject] :  (bool) plot average across all 12 subjets or not?
%   [saveFig]         :  (bool) save figures or not?
%
% Example:
%  makeFigure1('exampleSubject', 1, 'plotMeanSubject', false, 'saveFig', true)

%% 0. Set up paths and define parameters

p = inputParser;
p.KeepUnmatched = true;
p.addParameter('exampleSubject', 12);
p.addParameter('plotMeanSubject', false, @islogical)
p.addParameter('saveFig', false, @islogical);
p.parse(varargin{:});

% Rename variables
exampleSubject      = p.Results.exampleSubject;
plotMeanSubject     = p.Results.plotMeanSubject;
saveFig             = p.Results.saveFig;

% Define subjects
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

% What visual area to use?
area            = 'V123'; % Choose between 'V1', 'V2', 'V3', or 'V123'
eccenLimitDeg   = [2 6]; % what is the eccentricity limit (deg) for the template, supposingly matching the stimulus aperture. 
                       % (Can be a single int x, to get [0 x] or a vector limiting between [x,y])

% Define colormap and contour lines
contourmapPercentile   = 93.6; % draw contour line at what fraction of the colormap?  top 15 channels: 90.4, or for top 10 channels: 93.6, 
                                % or for use any integer under 10 to get contour lines at equal percentiles of data
colormapPercentile     = 97.5; % percentile of data to use for max/min limits of colorbar

% Number of iterations for the random coherence prediction of the forward model
n        	= 10;        % number of timepoints (ms)
nrEpochs    = 1000;      % number of epochs
theta       = 0;         % von mises mean, equal for three distributions
kappa.coh   = 10*pi;     % very narrow von Mises
kappa.incoh = 0;         % very broad (uniform) von Mises
kappa.mix   = 0.27*pi;        % in between width size von Mises

% Define vector that can truncate number of sensors 
keep_sensors = logical([ones(157,1); zeros(192-157,1)]); % TODO: Figure out a more generic way to define keep_sensors

% Loop over subjects
if plotMeanSubject
    subjectsToPlot = 1:length(subject);
else
    subjectsToPlot = exampleSubject;
end

for s = subjectsToPlot
    
    d = dir(fullfile(bsDB, projectName, 'data', subject{s}, 'R*'));
    bsData = fullfile(d(1).folder, d(1).name);    
    bsAnat = fullfile(bsDB, projectName, 'anat', subject{s});
    
    %% 1. Load relevant matrices
    
    G_constrained = getGainMatrix(bsData, keep_sensors);

    % Get V1 template limited to 11 degrees eccentricity
    template = getTemplate(bsAnat, area, eccenLimitDeg);

    % Simulate coherent, in between or mixture, adn incoherent source time 
    % series and compute predictions from forward model (w)
    tmp = getForwardModelPredictions(G_constrained, template.(area), [], n, nrEpochs, theta, kappa);
    
    % Compute amplitude across time
    amps.c = abs(fft(tmp.c,[],2));
    amps.i = abs(fft(tmp.i,[],2));
    amps.m = abs(fft(tmp.m,[],2));
    
    % Compute mean weights across epochs at input frequency
    w.V1c(s,:) = mean(amps.c(:,2,:),3);
    w.V1i(s,:) = mean(amps.i(:,2,:),3);
    w.V1m(s,:) = mean(amps.m(:,2,:),3);
    
end
    %% Visualize predictions

for exampleSubject = subjectsToPlot
    dataToPlot   = cat(1,w.V1c(exampleSubject,:), w.V1i(exampleSubject,:), w.V1m(exampleSubject,:));
    colorMarkers = {'r','b', 'g'};
    sub_ttl      = {sprintf('Uniform phase S%d', exampleSubject), ...
                    sprintf('Random phase S%d', exampleSubject),...
                    sprintf('Mixed phase S%d', exampleSubject)};                
    fig_ttl      = {sprintf('Figure1_model_predictions_mixture_%s_%d-%d', area, eccenLimitDeg(1),eccenLimitDeg(2)), sprintf('Figure1_Uniform_and_Random_Compared_mixture_%s_%d-%d', area, eccenLimitDeg(1),eccenLimitDeg(2))};
    markerType   = '.';

    dataDir       = fullfile(fmsRootPath,'data', subject{exampleSubject}); % Where to save images?
    figureDir       = fullfile(fmsRootPath,'figures', subject{exampleSubject}); % Where to save images?

    % Make figure and data dir for subject, if non-existing             
    if ~exist(figureDir,'dir'); mkdir(figureDir); end
    if ~exist(dataDir,'dir'); mkdir(dataDir); end

    sensorsWithinContours = visualizeSensormaps(dataToPlot, colormapPercentile, contourmapPercentile, colorMarkers, markerType, fig_ttl, sub_ttl, saveFig, figureDir);

    % Make them logicals so we can use them later as indices
%     sensorsWithinContours = logical(sensorsWithinContours);

    % Save sensors of interest falling within the contour lines
    if saveFig; save(fullfile(dataDir, sprintf('%s_prediction_%s_%d-%d.mat', subject{exampleSubject}, area, eccenLimitDeg(1),eccenLimitDeg(2))), 'dataToPlot'); end

end

%% Take mean across subjects and plot if requested
if plotMeanSubject
    
    w.V1c_mn = mean(w.V1c,1);
    w.V1i_mn = mean(w.V1i,1);
    w.V1m_mn = mean(w.V1m,1);
    
    dataToPlot = [w.V1c_mn; w.V1i_mn; w.V1m_mn];
    
    fig_ttl    = {sprintf('Figure1_model_predictions_mixture_%s_%d-%d', area, eccenLimitDeg(1),eccenLimitDeg(2)), sprintf('Figure1_Uniform_and_Random_Compared_mixture_%s_%d-%d', area, eccenLimitDeg(1),eccenLimitDeg(2))};
    sub_ttl    = {sprintf('Uniform phase Average N = %d', length(subject)), ...
                  sprintf('Random phase Average N = %d', length(subject)), ...
                  sprintf('Mixed phase Average N = %d', length(subject))};
    markerType = '.';
    
    figureDir       = fullfile(fmsRootPath,'figures', 'average'); % Where to save images?
    dataDir       = fullfile(fmsRootPath,'data', 'average'); % Where to save images?

    if ~exist(figureDir,'dir'); mkdir(figureDir); end
    if ~exist(dataDir,'dir'); mkdir(dataDir); end

    % Plot data and save channels that are located inside the contour lines
    sensorsWithinContours = visualizeSensormaps(dataToPlot, colormapPercentile, contourmapPercentile, colorMarkers, markerType, fig_ttl, sub_ttl, saveFig, figureDir);

    % Make them logicals so we can use them later as indices
%     sensorsWithinContours = logical(sensorsWithinContours);

    % Save sensors of interest falling within the contour lines
    if saveFig; save(fullfile(dataDir, sprintf('Average_prediction_%s_%d-%d.mat',area,eccenLimitDeg(1),eccenLimitDeg(2))), 'dataToPlot'); end       
            
end
