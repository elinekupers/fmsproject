function makeFigure6(varargin)

% This is a function to make Figure 6 from the manuscript using a forward
% model to predict coherent and incoherent neural sources to MEG responses, WITHOUT CANCELLATION.

% This figure shows the MEG forward model based on coherent and incoherent
% predictions coming from vertices located in V1 when using an ABSOLUTE gain matrix.

% To runs this script, you need:
% (1) Access to the SSMEG folder in the brainstorm data base
% (2) MEG_utils and Fieldtrip toolbox added to the paths. For example:
%     tbUse('ForwardModelSynchrony');
%        or to only add the MEG_utils toolbox:
%     addpath(genpath('~/matlab/git/toolboxes/meg_utils'))
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
p.addParameter('highResSurf', false, @islogical); % Use high or low resolution surface and plot ratio with vs without cancellation
p.parse(varargin{:});

% Rename variables
subjectsToPlot      = p.Results.subjectsToPlot;
plotMeanSubject     = p.Results.plotMeanSubject;
saveFig             = p.Results.saveFig;
highResSurf         = p.Results.highResSurf;

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
    'wlsubj070'}; %  S12 - From exp: Full only% Which subjects to show as example?
if nargin < 1; exampleSubject  = 12; end % Which example subject to show if not defined

% Path to brainstorm database and project name
bsDB            = '/Volumes/server/Projects/MEG/brainstorm_db/';
projectName     = 'SSMEG';


% What visual area to use?
area            = 'V123'; % Choose between 'V1', 'V2', 'V3', or 'V123'
eccenLimitDeg   = [.18 11]; % what is the eccentricity limit (deg) for the template, supposingly matching the stimulus aperture.
% (Can be a single int x, to get [0 x] or a vector [x,y] limiting eccentricity to larger/equal to x and smaller/equal to y)

% Define colormap and contour lines
contourmapPercentile   = 90.4; % draw contour line at what fraction of the colormap?  top 15 channels: 90.4, or for top 10 channels: 93.6,
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
end


for s = subjectsToPlot
    
    % Get subject data and anatomy files
    d = dir(fullfile(bsDB, projectName, 'data', subject{s}, 'R*'));
    bsData = fullfile(d(1).folder, d(1).name);
    
    if highResSurf
        bsAnat = fullfile(bsDB, projectName, 'anat', subject{s}, 'highres');
    else
        bsAnat = fullfile(bsDB, projectName, 'anat', subject{s});
    end
    
    %% 1. Load relevant matrices
    
    % Simulate coherent and incoherent source time series and compute
    % predictions from forward model (w)
    G_constrained = getGainMatrix(bsData, keep_sensors, highResSurf);
    
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
    w.woC.V1c(s,:) = mean(amps.woC.c(:,2,:),3);
    w.woC.V1i(s,:) = mean(amps.woC.i(:,2,:),3);
    
    w.wC.V1c(s,:) = mean(amps.wC.c(:,2,:),3);
    w.wC.V1i(s,:) = mean(amps.wC.i(:,2,:),3);
    
end

%% 3. Visualize predictions from forward model for requested individual subject

% Define plotting data and labels
dataToPlot   = cat(1, w.woC.V1c(exampleSubject,:), w.woC.V1i(exampleSubject,:));
colorMarkers = {'r','b'};
fig_ttl      = {'Figure6_V1_model_predictions-No_cancellation', ...
    'Figure6_Sl_and_Broadband_Compared-No_cancellation'};
sub_ttl      = {sprintf('No cancellation: Coherent phase S%d',exampleSubject), ...
    sprintf('No cancellation: Incoherent phase S%d',exampleSubject)};
markerType   = '.';

dataDir       = fullfile(fmsRootPath,'data', subject{s}); % Where to save vector of sensors that fall within contours?
figureDir       = fullfile(fmsRootPath,'figures', subject{s}); % Where to save images?

% Make figure and data dir for subject, if non-existing
if ~exist(figureDir,'dir'); mkdir(figureDir); end
if ~exist(dataDir,'dir'); mkdir(dataDir); end

% Plot it!
visualizeSensormaps(dataToPlot, colormapPercentile, contourmapPercentile, colorMarkers, markerType, fig_ttl, sub_ttl, saveFig, figureDir);

%% Visualize prediction for across subjects
if plotMeanSubject
    
    % Redefine figure dir
    figureDir    = fullfile(fmsRootPath, 'figures', 'average'); % Where to save images?
    
    % Take the average across subjects
    w.woC.V1c_mn     = mean(w.woC.V1c,1);
    w.woC.V1i_mn     = mean(w.woC.V1i,1);
    
    % Define plotting data and labels
    dataToPlot   = cat(1, w.woC.V1c_mn, w.woC.V1i_mn);
    colorMarkers = {'r','b'};
    fig_ttl      = {sprintf('Figure6_V123_model_predictions-No_cancellation_%s_%1.2f-%d_contour%03d_highResFlag%d', area, eccenLimitDeg(1),eccenLimitDeg(2), contourmapPercentile,highResSurf), ...
        sprintf('Figure6_Sl_and_Broadband_Compared-No_cancellation_%s_%1.2f-%d_contour%03d_highResFlag%d', area, eccenLimitDeg(1),eccenLimitDeg(2), contourmapPercentile,highResSurf)};
    sub_ttl      = {sprintf('No cancellation: Coherent phase Average N=%d', length(subject)), ...
        sprintf('No cancellation: Incoherent phase Average N=%d', length(subject))};
    markerType   = '.';
    
    % Plot it!
    visualizeSensormaps(dataToPlot, colormapPercentile, contourmapPercentile, colorMarkers, markerType, fig_ttl, sub_ttl, saveFig, figureDir);
end

%% Plot ratio of with vs without cancellation

% Take the average across subjects
w.wC.V1c_mn     = mean(w.wC.V1c,1);
w.wC.V1i_mn     = mean(w.wC.V1i,1);

% Take ratio
ratioCoh = w.woC.V1c_mn ./ w.wC.V1c_mn;
ratioInCoh = w.woC.V1i_mn ./ w.wC.V1i_mn;

clims = [-1 1].*10^-2;

figure;
subplot(311);
megPlotMap(w.wC.V1c_mn, 0.1*clims, [], 'bipolar', 'Coh: With cancellation');

subplot(312);
megPlotMap(w.woC.V1c_mn, 0.1*clims, [], 'bipolar', 'Coh: Without cancellation');

subplot(313);
megPlotMap(ratioCoh, [-1 1], [], 'bipolar', 'Coh: Ratio with / without');

figure;
subplot(311);
megPlotMap(w.wC.V1i_mn, 0.1*clims, [], 'bipolar', 'InCoh: With cancellation');

subplot(312);
megPlotMap(w.woC.V1i_mn, clims, [], 'bipolar', 'InCoh: Without cancellation');

subplot(313);
megPlotMap(ratioInCoh, [-1 1], [], 'bipolar', 'InCoh: Ratio with / without');

% Or just the two ratio's:
fig_ttl = {sprintf('Fig6_ratioWithVsWithoutCancellation_%s_%1.2f-%d_contour%d_highResFlag%d', area, eccenLimitDeg(1),eccenLimitDeg(2), contourmapPercentile,highResSurf), ...
    sprintf('Fig6_ratioWithVsWithoutCancellation_Overlap_%s_%1.2f-%d_contour%d_highResFlag%d', area, eccenLimitDeg(1),eccenLimitDeg(2), contourmapPercentile,highResSurf)};
sub_ttl = {'Coh: Ratio with / without', 'InCoh: Ratio with / without'};
visualizeSensormaps([ratioCoh; ratioInCoh], 100, [], colorMarkers, markerType, fig_ttl, sub_ttl, saveFig, figureDir);


return