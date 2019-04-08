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
area            = 'V123'; % Choose between 'V1', 'V2', 'V3', or 'V123'
eccenLimitDeg   = [2 6]; % what is the eccentricity limit (deg) for the template, supposingly matching the stimulus aperture.
                         % (Can be a single int x, to get [0 x] or a vector limiting between [x,y])

% Number of iterations for the random coherence prediction of the forward
% model
n               = 10;         % number of timepoints (ms)
nrEpochs        = 1000;       % number of epochs
theta           = 0;          % von mises mean of all three distributions
kappa.coh       = 10*pi;      % kappa width, coherent signal
kappa.incoh     = 0;          % kappa width, incoherent signal
allMixedKappas  = pi.*logspace(log10(.1),log10(2),10);

% Define vector that can truncate number of sensors
keep_sensors = logical([ones(157,1); zeros(192-157,1)]); % Note: Figure out a more generic way to define keep_sensors

% Plotting variables
nrMixedKappas = length(allMixedKappas);
labels = cellstr(sprintfc('Kappa = %1.2f *pi', allMixedKappas./pi));
clims =  10^-3.*[-1 1];
nrows = 4;
ncols = ceil(nrMixedKappas/nrows)+1;

% Two ways of plotting contours
% (1) one contour line at the percentile of data (say 90.4/100).
% (2) number of contour lines, dividing data into equal groups (use one number under 10)
%   for example, if contourPercentile=3, you draw 3 lines at the 25, 50 and 75th percentile
% contourPercentile = 90.4;
contourPercentile = 3;

% line with for contour lines
lw = 2;

for s = subjectToPlot
    
    d = dir(fullfile(bsDB, projectName, 'data', subject{s}, 'R*'));
    bsData = fullfile(d(1).folder, d(1).name);
    bsAnat = fullfile(bsDB, projectName, 'anat', subject{s});
    
    figure(1); set(1, 'Color', 'w', 'Position', [1, 1, 1680, 999]); clf; hold all;
    
    for k = 1:length(allMixedKappas)
        
        kappa.mix   = allMixedKappas(k);
        
        %% 1. Load relevant matrices
        
        % Get Gain matrix
        G_constrained = getGainMatrix(bsData, keep_sensors);
        
        % Get V1 template limited to 11 degrees eccentricity
        template = getTemplate(bsAnat, area, eccenLimitDeg);
        
        % Simulate coherent, in between or mixture, adn incoherent source time
        % series and compute predictions from forward model (w)
        tmp = getForwardModelPredictions(G_constrained, template.([area '_stimEccen']), [], n, nrEpochs, theta, kappa);
       
        
        % Compute amplitude across time
        amps.c = abs(fft(tmp.c,[],2));
        amps.i = abs(fft(tmp.i,[],2));
        amps.m = abs(fft(tmp.m,[],2));
        
        % Compute mean weights across epochs at input frequency
        w.V1c(s,k,:) = mean(amps.c(:,2,:),3);
        w.V1i(s,k,:) = mean(amps.i(:,2,:),3);
        w.V1m(s,k,:) = mean(amps.m(:,2,:),3);
        
    end
    
    % Visualize most incoherent signals
    subplot(nrows, ncols,1);
    dataToPlot = squeeze(w.V1i(s,1,:));
    if contourPercentile>10; contourLines = [1 1]*prctile(dataToPlot, contourPercentile); else contourLines = contourPercentile; end
    megPlotMap(dataToPlot, clims, [], 'bipolar', 'Kappa = 0', [],[], 'isolines', contourLines);
    h = findobj(gca,'Type','contour');
    h.LineWidth = lw;
    
    % Visualize most coherent signals
    subplot(nrows, ncols,nrMixedKappas+2);
    dataToPlot = squeeze(w.V1c(s,1,:));
    if contourPercentile>10; contourLines = [1 1]*prctile(dataToPlot, contourPercentile); else contourLines = contourPercentile; end
    megPlotMap(dataToPlot, clims, [], 'bipolar', 'Kappa = 10*pi', [],[], 'isolines', contourLines);
    h = findobj(gca,'Type','contour');
    h.LineWidth = lw;
    
    % Visualize all mixtures
    for k = 1:nrMixedKappas
        subplot(nrows, ncols,k+1)
        dataToPlot = squeeze(w.V1m(s,k,:));
        if contourPercentile>10; contourLines = [1 1]*prctile(dataToPlot, contourPercentile); else contourLines = contourPercentile; end
        megPlotMap(dataToPlot, clims, [], 'bipolar', labels(k), [],[], 'isolines', contourLines);
        h = findobj(gca,'Type','contour');
        h.LineWidth = lw;
    end
    
    dataDir       = fullfile(fmsRootPath,'data', subject{s}); % Where to save images?
    figureDir       = fullfile(fmsRootPath,'figures', subject{s}); % Where to save images?
    
    % Make figure and data dir for subject, if non-existing
    if ~exist(figureDir,'dir'); mkdir(figureDir); end
    if ~exist(dataDir,'dir'); mkdir(dataDir); end
    
    save(fullfile(dataDir, sprintf('%s_mixturePredictions_contour.mat', subject{s})), 'w');
    figurewrite(fullfile(figureDir, sprintf('%s_mixturePredictions_%s_%d-%d_%2.1f', subject{s}, area, eccenLimitDeg(1), eccenLimitDeg(2), contourPercentile)),[],[1 300],'.',1);
    
end


%% Take mean across subjects and plot if requested
if plotMeanSubject
    
    w.V1c_mn = mean(w.V1c,1);
    w.V1i_mn = mean(w.V1i,1);
    w.V1m_mn = mean(w.V1m,1);
    
    figure(1); set(gcf, 'Color', 'w', 'Position', [1, 1, 1680, 999]); clf; hold all;
    subplot(nrows, ncols,1);
    dataToPlot = w.V1i_mn(1,:);
    if contourPercentile>10; contourLines = [1 1]*prctile(dataToPlot, contourPercentile); else contourLines = contourPercentile; end
    megPlotMap(dataToPlot, clims, [], 'bipolar', 'Kappa = 0', [],[], 'isolines', contourLines);
    h = findobj(gca,'Type','contour');
    h.LineWidth = lw;
    
    subplot(nrows, ncols,nrMixedKappas+2);
    dataToPlot = w.V1c_mn(1,:);
    if contourPercentile>10; contourLines = [1 1]*prctile(dataToPlot, contourPercentile); else contourLines = contourPercentile; end
    megPlotMap(dataToPlot, clims, [], 'bipolar', 'Kappa = 10*pi', [],[], 'isolines', contourLines);
    h = findobj(gca,'Type','contour');
    h.LineWidth = lw;
    
    for k = 1:nrMixedKappas
        subplot(nrows, ncols,k+1);
        dataToPlot = w.V1m(k,:);
        if contourPercentile>10; contourLines = [1 1]*prctile(dataToPlot, contourPercentile); else contourLines = contourPercentile; end
        megPlotMap(dataToPlot, clims, [], 'bipolar', labels(k), [],[], 'isolines', contourLines);
        h = findobj(gca,'Type','contour');
        h.LineWidth = lw;
    end
    
    figureDir       = fullfile(fmsRootPath,'figures', 'average'); % Where to save images?
    dataDir       = fullfile(fmsRootPath,'data', 'average'); % Where to save images?
    
    if ~exist(figureDir,'dir'); mkdir(figureDir); end
    if ~exist(dataDir,'dir'); mkdir(dataDir); end
    
    % Plot data and save data
    save(fullfile(dataDir, 'mixturePredictions_averge.mat'), 'w');
    figurewrite(fullfile(figureDir, sprintf('mixturePredictions_%s_%d-%d_%2.1f', area, eccenLimitDeg(1), eccenLimitDeg(2),contourPercentile)),[],[1 300],'.',1);
    
end
