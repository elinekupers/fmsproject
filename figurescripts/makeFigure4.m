function makeFigure4()

% This is a function to make Figure 4 from the manuscript using a forward
% model to predict coherent and incoherent neural sources to MEG responses, WITHOUT CANCELLATION.

% This figure shows the MEG forward model based on coherent and incoherent
% predictions coming from vertices located in V1 when using an ABSOLUTE gain matrix.

% To runs this script, you need:
% (1) Access to the SSMEG folder in the brainstorm data base
% (2) MEG_utils and Fieldtrip toolbox added to the paths. For example:
%     tbUse('ForwardModelSynchrony');
%        or to only add the MEG_utils toolbox:
%     addpath(genpath('~/matlab/git/toolboxes/meg_utils'))

%% 0. Set up paths and define parameters

% Set up paths
figureDir       = fullfile(fmsRootPath, 'figures', subject{exampleSubject}); % Where to save images?
saveFigures     = true;         % Save figures in the figure folder?
plotMeanSubject = false;        % Plot average subject?
bsDB            = '/Volumes/server/Projects/MEG/brainstorm_db/'; % Path to brainstorm database

% Define project name, subject and data/anatomy folders
projectName = 'SSMEG';

% Which subjects to average?
%   Full  only: 'wlsubj048', 'wlsubj046','wlsubj039','wlsubj059', 'wlsubj067'
%   Full, Left, Right: 'wlsubj002','wlsubj004','wlsubj005','wlsubj006','wlsubj010','wlsubj011'
subject = {'wlsubj002','wlsubj004','wlsubj005','wlsubj006','wlsubj010','wlsubj011'};

% Which subjects to show as example?
exampleSubject  = 1;

% What visual area?
area    = 'V1'; % Choose from 'V1', or 'all' (V1-V3)

% What's the plotting range for individual example and average across
% subjects?
contourmapPercentile   = 93.6; % draw contour line at what fraction of the colormap?
colormapPercentile     = 97.5; % percentile of data to use for max/min limits of colorbar

% Number of iterations for the random coherence prediction of the forward
% model
n        = 10;     % number of timepoints (ms)
nrEpochs = 100;        % number of epochs

% Define vector that can truncate number of sensors 
keep_sensors = logical([ones(157,1); zeros(192-157,1)]); % Note: Figure out a more generic way to define keep_sensors

for s = 1:length(subject)
    
    d = dir(fullfile(bsDB, projectName, 'data', subject{s}, 'R*'));
    dataDir = fullfile(d(1).folder, d(1).name);    
    anatDir = fullfile(bsDB, projectName, 'anat', subject{s});
    
    %% 1. Load relevant matrices
    
    G_constrained = getGainMatrix(dataDir, keep_sensors);

    % Get V1 template limited to 11 degrees eccentricity
    template = getTemplate(anatDir, area, 11);

    % Simulate coherent and incoherent source time series and compute
    % predictions from forward model (w)

    % NOTE: Take absolute values of G_contrained - no cancellation possible
    tmp = getForwardModelPredictions(abs(G_constrained), template.V1StimEccen, [], n, nrEpochs);
   
    % Take mean amplitude across epochs
    amps.c = abs(fft(tmp.c,[],2));
    amps.i = abs(fft(tmp.i,[],2));
    
    w.V1c(s,:) = mean(amps.c(:,2,:),3);
    w.V1i(s,:) = mean(amps.i(:,2,:),3);
    
end

%% 3. Visualize predictions from forward model for requested individual subject

% Define plotting data and labels
dataToPlot   = cat(1, w.V1c(exampleSubject,:), w.V1i(exampleSubject,:));
colorMarkers = {'r','b'};
fig_ttl      = {'Figure4_V1_model_predictions-No_cancellation', ...
                'Figure4_Sl_and_Broadband_Compared-No_cancellation'};
sub_ttl      = {'Coherent phase S1', ...
                'Incoherent phase S2'};
markerType   = '.';

% Plot it!
visualizeSensormaps(dataToPlot, colormapPercentile, contourmapPercentile, colorMarkers, markerType, fig_ttl, sub_ttl, saveFigures, figureDir);

%% Visualize prediction for across subjects
if plotMeanSubject
    
    % Redefine figure dir
    figureDir    = fullfile(fmsRootPath, 'figures'); % Where to save images?

    % Take the average across subjects
    w.V1c_mn     = mean(w.V1c,1);
    w.V1i_mn     = mean(w.V1i,1);

    % Define plotting data and labels
    dataToPlot   = cat(1, w.V1c_mn, w.V1i_mn);
    colorMarkers = {'r','b'};
    fig_ttl      = {'Figure4_V1_model_predictions-No_cancellation', ...
                    'Figure4_Sl_and_Broadband_Compared-No_cancellation'};
    sub_ttl      = {sprintf('No cancellation: Coherent phase Average N=%d', length(subject)), ...
                    sprintf('No cancellation: Incoherent phase Average N=%d', length(subject))};
    markerType   = '.';

    % Plot it!
    visualizeSensormaps(dataToPlot, colormapPercentile, contourmapPercentile, colorMarkers, markerType, fig_ttl, sub_ttl, saveFigures, figureDir);

end
