function makeSupplementaryFigure9(varargin)
% This is a function to make model predictions for different mixtures of
% synchrony levels in early visual cortical sources as shown in
% Supplementary Figure 9 from the manuscript:
%
%   A visual encoding model links magnetoencephalography
%   signals to neural synchrony in human cortex.
%       by Kupers, Benson, Winawer (YEAR) JOURNAL.
%
% To runs this script, you need:
% (1) Access to the SSMEG folder in the brainstorm data base
% (2) MEG_utils and Fieldtrip toolbox added to the paths. For example:
%     tbUse('ForwardModelSynchrony');
%        or to only add the MEG_utils toolbox:
%     addpath(genpath('~/matlab/git/toolboxes/meg_utils'))
% (3) Run the s_visualAreasFS2BS script from this repository
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
%                                     'benson17', or 'wang15atlas' (note:
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
%   [signedColorbar]        :  (bool) true/false plot signed colormap or only
%                                     positive values.
%
% Example 1: Plot first subject
%  makeSupplementaryFigure9('subjectsToPlot', 1, 'plotMeanSubject', false, 'saveFig', true)
% Example 2: Plot example subject in manuscript (S12)
%  makeSupplementaryFigure9('subjectsToPlot', 12, 'plotMeanSubject', false, 'saveFig', true)
% Example 3: Plot all subjects and group average
%  makeSupplementaryFigure9('subjectsToPlot', 1:12, 'plotMeanSubject', true, 'saveFig', true)
%

p = inputParser;
p.KeepUnmatched = true;
p.addParameter('subjectsToPlot', 12);
p.addParameter('plotMeanSubject', true, @islogical)
p.addParameter('saveFig', true, @islogical);
p.addParameter('headmodelType', 'OS', @(x) any(validatestring(x,{'OS', 'BEM'})));
p.addParameter('highResSurf', false, @islogical);
p.addParameter('area', 'V123', @(x) any(validatestring(x,{'V1', 'V2', 'V3','V123', 'benson18atlas', 'wang15atlas'})));
p.addParameter('eccenLimitDeg', [0.18 11], @isnumeric);
p.addParameter('contourPercentile', 93.6, @isnumeric);
p.addParameter('signedColorbar', false, @islogical);
p.parse(varargin{:});

% Rename variables
subjectsToPlot        = p.Results.subjectsToPlot;
plotMeanSubject       = p.Results.plotMeanSubject;
saveFig               = p.Results.saveFig;
headmodelType         = p.Results.headmodelType;
highResSurf           = p.Results.highResSurf;
area                  = p.Results.area;
eccenLimitDeg         = p.Results.eccenLimitDeg;
contourPercentile     = p.Results.contourPercentile;
signedColorbar        = p.Results.signedColorbar;

%% 0. Define subjects and paths

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


% Load all subjects when plotting the mean
if plotMeanSubject
    subjectsToLoad = 1:length(subject);
else
    subjectsToLoad = subjectsToPlot;
end

% Path to brainstorm database and project name
bsDB            = '/Volumes/server/Projects/MEG/brainstorm_db/';
projectName     = 'SSMEG';

% Parameters for forward model predictions
n               = 10;         % number of timepoints (ms)
nrEpochs        = 1000;       % number of epochs
theta           = 0;          % von mises mean of all three distributions
kappa.syn       = 100*pi;     % kappa (von Mises width), coherent signal
kappa.asyn      = 0;          % kappa (von Mises width), incoherent signal
allMixedKappas  = pi.*logspace(log10(.1),log10(2),10); % kappa (von Mises width), mixed  signals

% Define vector that can truncate number of sensors
keep_sensors = logical([ones(157,1); zeros(192-157,1)]); % Note: Figure out a more generic way to define keep_sensors

% Plotting variables
nrMixedKappas = length(allMixedKappas);
labels    = cellstr(sprintfc('Kappa = %1.2f *pi', allMixedKappas./pi));
cmapData  = bipolar(64);

% Check limits of color map/bar
if ~signedColorbar
    cmapData = cmapData(ceil(length(cmapData)/2):end,:);
end

% Subplot layout
nrows = 4;
ncols = ceil(nrMixedKappas/nrows)+1;

% line width for contour lines
lw = 2;

%% Loop over subjects to get headmodels
for s = subjectsToLoad
    
    d = dir(fullfile(bsDB, projectName, 'data', subject{s}, 'R*'));
    bsData = fullfile(d(1).folder, d(1).name);
    
    if highResSurf
        bsAnat = fullfile(bsDB, projectName, 'anat', subject{s}, 'highres');
    else
        bsAnat = fullfile(bsDB, projectName, 'anat', subject{s}, 'lowres');
    end
        
    %% Loop over kappa params to get predictions
    for k = 1:length(allMixedKappas)
        
        % Define current kappa for intermediate syncrhony level
        kappa.mix   = allMixedKappas(k);
        
        %% 1. Load relevant matrices
        
        % Get gain matrix
        G_constrained = getGainMatrix(bsData, keep_sensors, headmodelType, highResSurf);
        
        % Get V1 template limited to 11 degrees eccentricity
        template = getTemplate(bsAnat, area, eccenLimitDeg);
        
        % Simulate synchronous, mixture, and asynchonous source time
        % series and compute predictions from forward model (w)
        tmp = getForwardModelPredictions(G_constrained, template.([area '_StimEccen']), [], n, nrEpochs, theta, kappa);
        
        % Compute amplitude across time
        amps.c = abs(fft(tmp.c,[],2)); % coherent/synchonous ampl at input freq
        amps.i = abs(fft(tmp.i,[],2)); % incoherent/asynchonous ..
        amps.m = abs(fft(tmp.m,[],2)); % mixed ..
        
        % Compute mean weights across epochs at input frequency
        w.V1c(s,k,:) = mean(amps.c(:,2,:),3);
        w.V1i(s,k,:) = mean(amps.i(:,2,:),3);
        w.V1m(s,k,:) = mean(amps.m(:,2,:),3);
        
    end
end

% Compute the max of the average fully synchronous prediction to normalize
% all other predicted amplitudes to
synMax100 = max(squeeze(nanmean(w.V1c(:,1,:),1)));

%% Visualize predictions
for s = subjectsToPlot

    figure(1); clf; set(1, 'Color', 'w', 'Position', [1, 1, 1680, 999]);

    clims = [0 max(squeeze(w.V1c(s,1,:))'./synMax100)];
    
    % ASYNCHRONOUS (KAPPA=0)
    subplot(nrows, ncols,1);
    dataToPlot = squeeze(w.V1i(s,1,:))'./synMax100;
    if contourPercentile>10; contourLines = [1 1].*prctile(dataToPlot, contourPercentile); else, contourLines = contourPercentile; end
    [~,ch] = megPlotMap(dataToPlot, clims, [], cmapData, 'Kappa = 0', [],[], 'isolines', contourLines);
    h = findobj(gca,'Type','contour');
    h.LineWidth = lw;
    set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);
    
    % SYNCHRONOUS (KAPPA=100*pi)
    subplot(nrows, ncols,nrMixedKappas+2);
    dataToPlot = squeeze(w.V1c(s,1,:))'./synMax100;
    if contourPercentile>10; contourLines = [1 1].*prctile(dataToPlot, contourPercentile); else contourLines = contourPercentile; end
    [~,ch] = megPlotMap(dataToPlot, clims, [], cmapData, 'Kappa = 100*pi', [],[], 'isolines', contourLines);
    h = findobj(gca,'Type','contour');
    h.LineWidth = lw;
    set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);

    % ALL MIXTURES
    for k = 1:nrMixedKappas
        subplot(nrows, ncols,k+1)
        dataToPlot = squeeze(w.V1m(s,k,:))'./synMax100;
        if contourPercentile>10; contourLines = [1 1].*prctile(dataToPlot, contourPercentile); else contourLines = contourPercentile; end
        [~,ch] = megPlotMap(dataToPlot, clims, [], cmapData, labels(k), [],[], 'isolines', contourLines);
        h = findobj(gca,'Type','contour');
        h.LineWidth = lw;
        set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);
    end
    
    % Save figure if requested
    if saveFig
        figureDir     = fullfile(fmsRootPath,'figures', subject{s}); % Where to save images?
        % Make figure dir for subject, if non-existing
        if ~exist(figureDir,'dir'); mkdir(figureDir); end
        % save as png for memory purposes
        figurewrite(fullfile(figureDir, sprintf('SupplFig9_mixturePredictions_%s_%d-%d_%2.1f_%s_highResFlag%d_S%d_%d', area, eccenLimitDeg(1), eccenLimitDeg(2), contourPercentile, headmodelType, highResSurf, s, k)),0,[1 300],'.',1);
    end
    
end


%% Take mean across subjects and plot if requested
if plotMeanSubject
    
    % Get contour line colors
    colorContourLines = [linspace(57, 239, 12); linspace(83, 160, 12); linspace(164, 34, 12)];
    colorContourLines = colorContourLines'./255;
    
    % set colormap limits
    clims = [0 1];
    
    % Take average across subjects
    w.V1c_mn = squeeze(nanmean(w.V1c,1));
    w.V1i_mn = squeeze(nanmean(w.V1i,1));
    w.V1m_mn = squeeze(nanmean(w.V1m,1));
    
    fH1 = figure(1); clf; set(gcf, 'Color', 'w', 'Position', [1, 1, 1680, 999]); 
    fH2 = figure(2); clf; set(gcf, 'Color', 'w'); megPlotMap(zeros(1,157));
    
    %  Plot ASYNCHRONOUS/INCOHERENT (kappa = 0)
    figure(fH1); hold all;
    subplot(nrows, ncols,1);
    dataToPlot = w.V1i_mn(1,:)./synMax100;
    if contourPercentile>10; contourLines = [1 1].*prctile(dataToPlot, contourPercentile); else contourLines = contourPercentile; end
    [~,ch] = megPlotMap(dataToPlot, clims, [], cmapData, 'Kappa = 0', [],[], 'isolines', contourLines);
    h = findobj(gca,'Type','contour');
    h.LineWidth = lw;
    set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);

    figure(fH2); hold all;
    contourf(h.XData, h.YData, h.ZData, contourLines, 'LineColor',colorContourLines(1,:), 'Fill','off','LineWidth',2);
    colorbar off;

    % Plot SYNCHRONOUS/COHERENT (kappa = 100*pi)
    figure(fH1); 
    subplot(nrows, ncols,nrMixedKappas+2);
    dataToPlot = w.V1c_mn(1,:)./synMax100;
    if contourPercentile>10; contourLines = [1 1].*prctile(dataToPlot, contourPercentile); else contourLines = contourPercentile; end
    [~,ch] = megPlotMap(dataToPlot, clims, [], cmapData, 'Kappa = 100*pi', [],[], 'isolines', contourLines);
    h = findobj(gca,'Type','contour');
    h.LineWidth = lw;
    set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);
    
    % Plot contour lines into the single mesh
    figure(fH2); hold all;
    contourf(h.XData, h.YData, h.ZData, contourLines, 'LineColor',colorContourLines(12,:), 'Fill','off','LineWidth',2);
    colorbar off;
    
    % Plot meshes for MIXTURE Synchrony
    for k = 1:nrMixedKappas
        figure(fH1);
        subplot(nrows, ncols,k+1);
        dataToPlot = w.V1m_mn(k,:)./synMax100;
        if contourPercentile>10; contourLines = [1 1].*prctile(dataToPlot, contourPercentile); else contourLines = contourPercentile; end
        [~,ch] = megPlotMap(dataToPlot, clims, [], cmapData, labels(k), [], [], 'isolines', contourLines);
        h = findobj(gca,'Type','contour');
        h.LineWidth = lw;
        set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);
        
        % Plot contour lines into the single mesh
        figure(fH2); hold all; 
        contourf(h.XData, h.YData, h.ZData, contourLines, 'LineColor',colorContourLines(k+1,:), 'Fill','off','LineWidth',2);
        colorbar off;
    end
        
    % Save figure if requested
    if saveFig
        % Where to save images?
        figureDir  = fullfile(fmsRootPath,'figures', 'average');   
        if ~exist(figureDir,'dir'); mkdir(figureDir); end
        figure(fH1); % save as .png for memory purposes
        figurewrite(fullfile(figureDir, sprintf('SupplFig9_mixturePredictions_%s_%d-%d_%2.1f_%s_highResFlag%d_AVERAGE', area, eccenLimitDeg(1), eccenLimitDeg(2),contourPercentile, headmodelType, highResSurf)),0,[1 300],'.',1);
        figure(fH2);
        figurewrite(fullfile(figureDir, sprintf('SupplFig9_Contour_mixturePredictions_%s_%d-%d_%2.1f_%s_highResFlag%d_AVERAGE', area, eccenLimitDeg(1), eccenLimitDeg(2),contourPercentile, headmodelType, highResSurf)),[],0,'.',1);
    end
    
end
