function makeSupplementaryFigureXX_overlapDataWithPredictions()

% This is a function to make Supplementary Figure 1 from the manuscript about forward
% modeling coherent and incoherent neural sources to MEG responses.

% This figure shows the MEG forward model based on coherent and incoherent
% predictions coming from vertices located in V1, V2 and V3.

% To runs this script, you need:     
% (1) Access to the SSMEG folder in the brainstorm data base
% (2) MEG_utils and Fieldtrip toolbox added to the paths. For example:
%     tbUse('ForwardModelSynchrony');
%        or to only add the MEG_utils toolbox:
%     addpath(genpath('~/matlab/git/toolboxes/meg_utils'))


%% 0. Set up paths and define parameters

% Path to brainstorm database
bsDB            = '/Volumes/server/Projects/MEG/brainstorm_db/';
figureDir       = fullfile(fmsRootPath, 'figures'); % Where to save images?
saveFigures     = true;     % Save figures in the figure folder?

% Define project name, subject and data/anatomy folders
projectName    = 'SSMEG';

% Which subjects to average?
subject         = {'wl_subj002','wl_subj004','wl_subj005','wl_subj006','wl_subj010','wl_subj011'};
exampleSubject  = 2;

% What's the plotting range for individual example and average across
% subjects?
contourmapPercentile   = 90; % draw contour line at what fraction of the colormap?
colormapPercentile     = 97.5; % percentile of data to use for max/min limits of colorbar

% Number of iterations for the random coherence prediction of the forward
% model
n        = 10;         % number of timepoints (ms)
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
    template = getTemplate(anatDir, 'all', 11);

    % Simulate coherent and incoherent source time series and compute
    % predictions from forward model (w)
    tmp = getForwardModelPredictions(abs(G_constrained), template.V123StimEccen, [], n, nrEpochs);
   
    % Take mean amplitude across epochs
    amps.c = abs(fft(tmp.c,[],2));
    amps.i = abs(fft(tmp.i,[],2));
    
    w.V1c(s,:) = mean(amps.c(:,2,:),3);
    w.V1i(s,:) = mean(amps.i(:,2,:),3);
    
end

%% Take mean across subjects

w.V1c_mn = mean(w.V1c,1);
w.V1i_mn = mean(w.V1i,1);

%% Visualize predictions

dataAll      = cat(1,w.V1c(exampleSubject,:), w.V1i(exampleSubject,:), w.V1c_mn, w.V1i_mn);
colorMarkers = {'r','b', 'r', 'b'};
sub_ttl      = {sprintf('Coherent phase S%d', exampleSubject), ...
                sprintf('Incoherent phase S%d', exampleSubject), ...
                'Coherent phase Average S1-S6', ...
                'Incoherent phase Average S1-S6'};
fig_ttl      = {'SupplFigure1_V123_model_predictions', 'SupplFigure1_V123Coherent_and_Incoherent_Compared'};
markerType   = '.';

visualizeSensormaps(dataAll, colormapPercentile, contourmapPercentile, colorMarkers, markerType, fig_ttl, sub_ttl, saveFigures, figureDir);
