function makeFigure4(exampleSubject)

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

% Which subjects to average?
%   Full  only: 'wlsubj048', 'wlsubj046','wlsubj039','wlsubj059','wlsubj067',' wlsubj070'
%   Full, Left, Right: 'wlsubj002','wlsubj004','wlsubj005','wlsubj006','wlsubj010','wlsubj011'
subject = {'wlsubj002','wlsubj004','wlsubj005','wlsubj006','wlsubj010','wlsubj011','wlsubj048', 'wlsubj046','wlsubj039','wlsubj059', 'wlsubj067','wlsubj070'};

% Which subjects to show as example?
if nargin < 1; exampleSubject  = 12; end % Which example subject to show if not defined

% Path to brainstorm database and project name
bsDB            = '/Volumes/server/Projects/MEG/brainstorm_db/';
projectName     = 'SSMEG';
figureDir       = fullfile(fmsRootPath,'figures', subject{exampleSubject}); % Where to save images?
% dataDir         = fullfile(fmsRootPath,'data', subject{exampleSubject}); % Where to save images?
saveFigures     = true;     % Save figures in the figure folder?
plotMeanSubject = true;     % Plot average subject?
plotWithVsWithoutCancellation = true;

% What visual area?
area    = 'all'; % Choose from 'V1', or 'all' (V1-V3)

% What's the plotting range for individual example and average across
% subjects?
contourmapPercentile   = 93.6; %Choose 90.4 for top 15, or 93.6 for top 10; % draw contour line at what fraction of the colormap?
colormapPercentile     = 97.5; % percentile of data to use for max/min limits of colorbar

% Number of iterations for the random coherence prediction of the forward
% model
n        = 10;     % number of timepoints (ms)
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

    % NOTE: Take absolute values of G_contrained - no cancellation possible
    if strcmp(area, 'all')
        tmp.woC = getForwardModelPredictions(abs(G_constrained), template.V123StimEccen, [], n, nrEpochs);
        tmp.wC = getForwardModelPredictions(G_constrained, template.V123StimEccen, [], n, nrEpochs);
    else
       tmp.woC = getForwardModelPredictions(abs(G_constrained), template.V1StimEccen, [], n, nrEpochs);
       tmp.wC = getForwardModelPredictions(G_constrained, template.V1StimEccen, [], n, nrEpochs);
    end
    
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
fig_ttl      = {'Figure4_V1_model_predictions-No_cancellation', ...
                'Figure4_Sl_and_Broadband_Compared-No_cancellation'};
sub_ttl      = {sprintf('No cancellation: Coherent phase S%d',exampleSubject), ...
                sprintf('No cancellation: Incoherent phase S%d',exampleSubject)};
markerType   = '.';

% Plot it!
visualizeSensormaps(dataToPlot, colormapPercentile, contourmapPercentile, colorMarkers, markerType, fig_ttl, sub_ttl, saveFigures, figureDir);

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
    fig_ttl      = {'Figure4_V1_model_predictions-No_cancellation', ...
                    'Figure4_Sl_and_Broadband_Compared-No_cancellation'};
    sub_ttl      = {sprintf('No cancellation: Coherent phase Average N=%d', length(subject)), ...
                    sprintf('No cancellation: Incoherent phase Average N=%d', length(subject))};
    markerType   = '.';

    % Plot it!
    visualizeSensormaps(dataToPlot, colormapPercentile, contourmapPercentile, colorMarkers, markerType, fig_ttl, sub_ttl, saveFigures, figureDir);

    
    %% Plot ratio of with vs without cancellation
    if plotWithVsWithoutCancellation

        % Take the average across subjects
        w.wC.V1c_mn     = mean(w.wC.V1c,1);
        w.wC.V1i_mn     = mean(w.wC.V1i,1);
                
        % Take ratio
        ratioCoh = w.wC.V1c_mn ./ w.woC.V1c_mn;
        ratioInCoh = w.wC.V1i_mn ./ w.woC.V1i_mn;

       clims = [-1 1].*10^-2;

       figure; 
       subplot(311); 
       megPlotMap(w.wC.V1c_mn, 0.1*clims, [], 'bipolar', 'Coh: With cancellation'); 
       
       subplot(312); 
       megPlotMap(w.woC.V1c_mn, 2*clims, [], 'bipolar', 'Coh: Without cancellation'); 
       
       subplot(313); 
       megPlotMap(ratioCoh, [-1 1], [], 'bipolar', 'Coh: Ratio with / without'); 
       
       figure; 
       subplot(311); 
       megPlotMap(w.wC.V1i_mn, 0.1*clims, [], 'bipolar', 'InCoh: With cancellation'); 
       
       subplot(312); 
       megPlotMap(w.woC.V1i_mn, 0.1*clims, [], 'bipolar', 'InCoh: Without cancellation'); 
       
       subplot(313); 
       megPlotMap(ratioInCoh, [-1 1], [], 'bipolar', 'InCoh: Ratio with / without'); 

        % Or just the two ratio's:
       fig_ttl = {'Fig4_ratioWithVsWithoutCancellation','Fig4_ratioWithVsWithoutCancellation_Overlap'};
       sub_ttl = {'Coh: Ratio with / without', 'InCoh: Ratio with / without'};
       visualizeSensormaps([ratioCoh; ratioInCoh], 100, contourmapPercentile, colorMarkers, markerType, fig_ttl, sub_ttl, saveFigures, figureDir);

    end

end