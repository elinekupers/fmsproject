function makeFigure1Mixtures(subjectToPlot)

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

if ~exist('subjectToPlot', 'var') || isempty(subjectToPlot)
    subjectToPlot = 1:12;
    plotMeanSubject = true;     % Plot average subject?
else
    plotMeanSubject = false;     % Plot average subject?
end

%% 0. Set up paths and define parameters

% Which subjects to average?
%   Full  only: 'wlsubj048', 'wlsubj046','wlsubj039','wlsubj059', 'wlsubj067'
%   Full, Left, Right: 'wlsubj002','wlsubj004','wlsubj005','wlsubj006','wlsubj010','wlsubj011'
subject         = {'wlsubj002','wlsubj004','wlsubj005','wlsubj006','wlsubj010','wlsubj011','wlsubj048', 'wlsubj046','wlsubj039','wlsubj059', 'wlsubj067', 'wlsubj070'};

% Path to brainstorm database and project name
bsDB            = '/Volumes/server/Projects/MEG/brainstorm_db/';
projectName     = 'SSMEG';
saveFigures     = true;     % Save figures in the figure folder?

% What visual area to use?
area            = 'all'; % Choose between 'V1' or 'all' (=V1-V3);

% Number of iterations for the random coherence prediction of the forward
% model
n           = 10;         % number of timepoints (ms)
nrEpochs    = 1000;       % number of epochs
theta       = 0;          % von mises mean of three distributions
kappa.coh   = 10*pi;
kappa.incoh = 0;
allMixedKappas   = [pi*0.125, pi.*linspace(0.25, 5, 20)];

% Define vector that can truncate number of sensors
keep_sensors = logical([ones(157,1); zeros(192-157,1)]); % Note: Figure out a more generic way to define keep_sensors

% Plotting variables
nrMixedKappas = length(allMixedKappas);
labels = cellstr(sprintfc('Kappa = %1.2f *pi', allMixedKappas./pi));
clims =  10^-3.*[-1 1];
nrows = 3;
ncols = ceil(nrMixedKappas/3)+1;

for s = subjectToPlot
    
     d = dir(fullfile(bsDB, projectName, 'data', subject{s}, 'R*'));
     bsData = fullfile(d(1).folder, d(1).name);
     bsAnat = fullfile(bsDB, projectName, 'anat', subject{s});
     
     figure(1); set(gcf, 'Color', 'w', 'Position', [1, 1, 1680, 999]); clf; hold all;
    
    for k = 1:length(allMixedKappas)
        
        kappa.mix   = allMixedKappas(k);
              
        %% 1. Load relevant matrices
        
        G_constrained = getGainMatrix(bsData, keep_sensors);
        
        % Get V1 template limited to 11 degrees eccentricity
        template = getTemplate(bsAnat, area, 11);
        
        % Simulate coherent, in between or mixture, adn incoherent source time
        % series and compute predictions from forward model (w)
        if strcmp(area, 'all')
            tmp = getForwardModelPredictions(G_constrained, template.V123StimEccen, [], n, nrEpochs, theta, kappa);
        else
            tmp = getForwardModelPredictions(G_constrained, template.V1StimEccen, [], n, nrEpochs, theta, kappa);
        end
        
        % Compute amplitude across time
        amps.c = abs(fft(tmp.c,[],2));
        amps.i = abs(fft(tmp.i,[],2));
        amps.m = abs(fft(tmp.m,[],2));
        
        % Compute mean weights across epochs at input frequency
        w.V1c(s,k,:) = mean(amps.c(:,2,:),3);
        w.V1i(s,k,:) = mean(amps.i(:,2,:),3);
        w.V1m(s,k,:) = mean(amps.m(:,2,:),3);
        
    end
    
    % Visualize predictions
    subplot(nrows, ncols,1);
    megPlotMap(squeeze(w.V1i(s,1,:)), clims, [], 'bipolar', 'Kappa = 0');
    
    subplot(nrows, ncols,nrMixedKappas+2);
    megPlotMap(squeeze(w.V1c(s,1,:)), clims, [], 'bipolar', 'Kappa = 10*pi');

    for k = 1:nrMixedKappas
        subplot(nrows, ncols,k+1)
        megPlotMap(squeeze(w.V1m(s,k,:)), clims, [], 'bipolar', labels(k));
    end
        
    dataDir       = fullfile(fmsRootPath,'data', subject{s}); % Where to save images?
    figureDir       = fullfile(fmsRootPath,'figures', subject{s}); % Where to save images?
    
    % Make figure and data dir for subject, if non-existing
    if ~exist(figureDir,'dir'); mkdir(figureDir); end
    if ~exist(dataDir,'dir'); mkdir(dataDir); end
        
    save(fullfile(dataDir, sprintf('%s_mixturePredictions.mat', subject{s})), 'w');
    figurewrite(fullfile(figureDir, sprintf('%s_mixturePredictions_%s', subject{s}, area)),[],0,'.',1);

end


%% Take mean across subjects and plot if requested
if plotMeanSubject
    
    w.V1c_mn = mean(w.V1c,1);
    w.V1i_mn = mean(w.V1i,1);
    w.V1m_mn = mean(w.V1m,1);
    
    figure(1); set(gcf, 'Color', 'w', 'Position', [1, 1, 1680, 999]); clf; hold all;
    subplot(nrows, ncols,1);
    megPlotMap(w.V1c_mn(1,:), clims, [], 'bipolar', 'Kappa = 0');
    
    subplot(nrows, ncols,nrMixedKappas+2);
    megPlotMap(w.V1i_mn(1,:), clims, [], 'bipolar', 'Kappa = 10*pi');

    for k = 1:nrMixedKappas
        subplot(nrows, ncols,k+1)
        megPlotMap(w.V1m(k,:), clims, [], 'bipolar', labels(k));
    end
        
    figureDir       = fullfile(fmsRootPath,'figures', 'average'); % Where to save images?
    dataDir       = fullfile(fmsRootPath,'data', 'average'); % Where to save images?
    
    if ~exist(figureDir,'dir'); mkdir(figureDir); end
    if ~exist(dataDir,'dir'); mkdir(dataDir); end
    
    % Plot data and save data
    save(fullfile(dataDir, 'mixturePredictions_averge.mat'), 'w');
    figurewrite(fullfile(figureDir, sprintf('mixturePredictions_%s', area)),[],0,'.',1);

end
