function makeFigure1(exampleSubject)

% This is a function to make Figure 1 from the manuscript about forward
% modeling coherent and incoherent neural sources to MEG responses.

% This figure shows the MEG forward model based on coherent and incoherent
% predictions coming from vertices located in V1.

% To runs this script, you need:     
% (1) Access to the SSMEG folder in the brainstorm data base
% (2) MEG_utils and Fieldtrip toolbox added to the paths. For example:
%     tbUse('ForwardModelSynchrony');
%        or to only add the MEG_utils toolbox:
%     addpath(genpath('~/matlab/git/toolboxes/meg_utils'))
% (3) Run the s_visualAreasFS2BS script from this repository


%% 0. Set up paths and define parameters

% Which subjects to average?
%   Full  only: 'wlsubj048', 'wlsubj046','wlsubj039','wlsubj059', 'wlsubj067'
%   Full, Left, Right: 'wlsubj002','wlsubj004','wlsubj005','wlsubj006','wlsubj010','wlsubj011'
subject         = {'wlsubj002','wlsubj004','wlsubj005','wlsubj006','wlsubj010','wlsubj011','wlsubj048', 'wlsubj046','wlsubj039','wlsubj059', 'wlsubj067', 'wlsubj070'};
if nargin < 1; exampleSubject  = 12; end

% Path to brainstorm database and project name
bsDB            = '/Volumes/server/Projects/MEG/brainstorm_db/';
projectName     = 'SSMEG';
figureDir       = fullfile(fmsRootPath,'figures', subject{exampleSubject}); % Where to save images?
dataDir         = fullfile(fmsRootPath,'data', subject{exampleSubject}); % Where to save images?
saveFigures     = false;     % Save figures in the figure folder?
plotMeanSubject = true;     % Plot average subject?

% What visual area to use?
area            = 'all'; % Choose between 'V1' or 'all' (=V1-V3);

% What's the plotting range for individual example and average across
% subjects?
contourmapPercentile   = 93.6;%93.6; % top 15: 90.4, or for top 10: 93.6; % draw contour line at what fraction of the colormap?
colormapPercentile     = 97.5; % percentile of data to use for max/min limits of colorbar

% Number of iterations for the random coherence prediction of the forward
% model
n        = 10;         % number of timepoints (ms)
nrEpochs = 1000;        % number of epochs

% Define vector that can truncate number of sensors 
keep_sensors = logical([ones(157,1); zeros(192-157,1)]); % Note: Figure out a more generic way to define keep_sensors


for s = 1:length(subject)
    
    d = dir(fullfile(bsDB, projectName, 'data', subject{s}, 'R*'));
    bsData = fullfile(d(1).folder, d(1).name);    
    bsAnat = fullfile(bsDB, projectName, 'anat', subject{s});
    
    %% 1. Load relevant matrices
    
    G_constrained = getGainMatrix(bsData, keep_sensors);

    % Get V1 template limited to 11 degrees eccentricity
    template = getTemplate(bsAnat, area, 11);

    % Simulate coherent and incoherent source time series and compute
    % predictions from forward model (w)
    if strcmp(area, 'all')
        tmp = getForwardModelPredictions(G_constrained, template.V123StimEccen, [], n, nrEpochs);
    else
       tmp = getForwardModelPredictions(G_constrained, template.V1StimEccen, [], n, nrEpochs);
    end
    
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
for exampleSubject = 1:12
    dataToPlot   = cat(1,w.V1c(exampleSubject,:), w.V1i(exampleSubject,:), w.V1m(exampleSubject,:));
    colorMarkers = {'r','b', 'g'};
    sub_ttl      = {sprintf('Uniform phase S%d', exampleSubject), ...
                    sprintf('Random phase S%d', exampleSubject),...
                    sprintf('Mixed phase S%d', exampleSubject)};                
    fig_ttl      = {sprintf('Figure1_model_predictions_mixture_%s', area), sprintf('Figure1_Uniform_and_Random_Compared_mixture_%s', area)};
    markerType   = '.';

    % Make figure and data dir for subject, if non-existing             
    if ~exist(figureDir,'dir'); mkdir(figureDir); end
    if ~exist(dataDir,'dir'); mkdir(dataDir); end

    sensorsOfInterest = visualizeSensormaps(dataToPlot, colormapPercentile, contourmapPercentile, colorMarkers, markerType, fig_ttl, sub_ttl, saveFigures, figureDir);

    % Make them logicals so we can use them later as indices
%     sensorsOfInterest = logical(sensorsOfInterest);

    % Save sensors of interest falling within the contour lines
%     save(fullfile(dataDir, sprintf('%s_sensorsOfInterestFromPrediction.mat', subject{exampleSubject})), 'sensorsOfInterest');
%     save(fullfile(dataDir, sprintf('%s_prediction.mat', subject{exampleSubject})), 'dataToPlot');

end


%% Take mean across subjects and plot if requested
if plotMeanSubject
    
    w.V1c_mn = mean(w.V1c,1);
    w.V1i_mn = mean(w.V1i,1);
    w.V1m_mn = mean(w.V1m,1);
    
    dataToPlot = [w.V1c_mn; w.V1i_mn; w.V1m_mn];
    
    fig_ttl    = {sprintf('Figure1_model_predictions_mixture_%s', area), sprintf('Figure1_Uniform_and_Random_Compared_mixture_%s', area)};
    sub_ttl    = {sprintf('Uniform phase Average N = %d', length(subject)), ...
                  sprintf('Random phase Average N = %d', length(subject)),...
                  sprintf('Mixed phase Average N = %d', length(subject))};
    markerType = '.';
    
    figureDir       = fullfile(fmsRootPath,'figures', 'average'); % Where to save images?
    dataDir       = fullfile(fmsRootPath,'data', 'average'); % Where to save images?

    if ~exist(figureDir,'dir'); mkdir(figureDir); end
    if ~exist(dataDir,'dir'); mkdir(dataDir); end

    % Plot data and save channels that are located inside the contour lines
    sensorsOfInterest = visualizeSensormaps(dataToPlot, colormapPercentile, contourmapPercentile, colorMarkers, markerType, fig_ttl, sub_ttl, saveFigures, figureDir);

    % Make them logicals so we can use them later as indices
    sensorsOfInterest = logical(sensorsOfInterest);

    % Save sensors of interest falling within the contour lines
    save(fullfile(dataDir, sprintf('Average_sensorsOfInterestFromPrediction_%s', area)), 'sensorsOfInterest');
    save(fullfile(dataDir, sprintf('Average_prediction_%s',area)), 'dataToPlot');           
            
end
