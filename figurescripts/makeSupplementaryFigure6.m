function makeSupplementaryFigure6(varargin)
% This is a function to make Supplementary Figure 6 from the manuscript:
%   A visual encoding model links magnetoencephalography
%   signals to neural synchrony in human cortex.
%       by Kupers, Benson, Winawer (YEAR) JOURNAL.

% This figure shows model predictions for different stimulus eccentricities
% for asynchronous and synchronous sources in early visual cortical
% sources.
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
%                                     Choose from 'V1', 'V2', 'V3', 'V123',
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
%   [singleColorbarFlag]    :  (bool) Use a single colorbar per row,instead of 
%                                     per individual plot (default: true)
%
% Example 1:
%  makeSupplementaryFigure6('subjectsToPlot', 1, 'plotMeanSubject', false, 'saveFig', true)
% Example 2:
%  makeSupplementaryFigure6('subjectsToPlot', 12, 'plotMeanSubject', false, 'saveFig', true)
% Example 3:
%  makeSupplementaryFigure6('subjectsToPlot', 1:12, 'plotMeanSubject', true, 'saveFig', true)
%

p = inputParser;
p.KeepUnmatched = true;
p.addParameter('subjectsToPlot', 12);
p.addParameter('plotMeanSubject', true, @islogical)
p.addParameter('saveFig', true, @islogical);
p.addParameter('headmodelType', 'OS', @(x) any(validatestring(x,{'OS', 'BEM'})));
p.addParameter('highResSurf', false, @islogical);
p.addParameter('area', 'V123', @(x) any(validatestring(x,{'V1', 'V2', 'V3','V123', 'benson18atlas', 'wang15atlas'})));
p.addParameter('eccenLimitDeg', [[0.18,1]; [0.18,2]; [0.18,4]; [0.18,8]; [0.18,11]], @isnumeric);
p.addParameter('contourPercentile', 93.6, @isnumeric);
p.addParameter('signedColorbar', false, @islogical);
p.addParameter('singleColorbarFlag', true, @islogical);
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

% Get subject names and corresponding data session number
subject = getSubjectIDs;

% Load all subjects when plotting the mean
if plotMeanSubject
    subjectsToLoad = 1:length(subject);
else
    subjectsToLoad = subjectsToPlot;
end

% Path to brainstorm database and project name
bsDB            = '/Volumes/server/Projects/MEG/brainstorm_db/';
projectName     = 'SSMEG';

% Forward model time series specifics
n               = 10;         % number of timepoints (ms)
nrEpochs        = 1000;       % number of epochs
theta           = 0;          % von mises mean of all three distributions
kappa.syn       = 100*pi;      % kappa width, coherent signal
kappa.asyn      = 0;          % kappa width, incoherent signal
kappa.mix       = 0;

% Define vector that can truncate number of sensors
keep_sensors = logical([ones(157,1); zeros(192-157,1)]); % Note: Figure out a more generic way to define keep_sensors

% Plotting variables
nrEccenLimits = size(eccenLimitDeg,1);
for ii = 1:nrEccenLimits
    labels{ii}    = sprintf('Eccen = %1.2f-%1.2f', eccenLimitDeg(ii,:));
end

% Check limits of color map/bar
cmapData  = bipolar(64);
if ~signedColorbar
    cmapData = cmapData(ceil(length(cmapData)/2):end,:);
end

% Subplot layout
nrows = 2;
ncols = nrEccenLimits;

% line with for contour lines
lw = 2;
colorContourLines = jet(nrEccenLimits);

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
    for el = 1:nrEccenLimits
        
        thisEccen = eccenLimitDeg(el,:);
        
        %% 1. Load relevant matrices
        
        G_constrained = getGainMatrix(bsData, keep_sensors, headmodelType, highResSurf);
        
        % Get V1 template limited to 11 degrees eccentricity
        template = getTemplate(bsAnat, area, thisEccen);
        
        % Simulate coherent, in between or mixture, adn incoherent source time
        % series and compute predictions from forward model (w)
        tmp = getForwardModelPredictions(G_constrained, template.([area '_StimEccen']), [], n, nrEpochs, theta, kappa);
        
        % Compute amplitude across time
        amps.c = abs(fft(tmp.c,[],2));
        amps.i = abs(fft(tmp.i,[],2));
        
        % Compute mean weights across epochs at input frequency
        w.V1c(s,el,:) = mean(amps.c(:,2,:),3);
        w.V1i(s,el,:) = mean(amps.i(:,2,:),3);
        
        
    end
end

% Get max of synchronous full stim group average, to normalize predictions
maxSynFullStim = squeeze(max(nanmean(w.V1c(:,5,:),1)));


%% Visualize predictions
for s = subjectsToPlot
    
    fH1 = figure(1); set(fH1, 'Color', 'w', 'Position', [1, 1, 1680, 999]); clf; hold all;
    fH2 = figure(2); clf; set(gcf, 'Color', 'w'); hold all; megPlotMap(zeros(1,157)); drawnow; title('Synchronous sources')
    fH3 = figure(3); clf; set(gcf, 'Color', 'w'); hold all; megPlotMap(zeros(1,157)); drawnow; title('Asynchronous sources')
    
    colormapLimits = [0 max(squeeze(w.V1c(s,nrEccenLimits,:))./maxSynFullStim)];
    
    for el = 1:nrEccenLimits
                
        % SYNCHRONOUS (KAPPA=100*pi)
        figure(fH1); hold all;
        subplot(nrows, ncols, el);
        dataToPlot = squeeze(w.V1c(s,el,:))./maxSynFullStim;

        % Check at what point to draw contour lines
        if ~isempty(contourPercentile)
            if contourPercentile > 10 % Plot at given data percentile
                contourLines = [1 1]*prctile(dataToPlot, contourPercentile);
            else % Plot nr of contours at equo-distance percentiles
                contourLines = contourPercentile;
            end
        else, contourLines = []; end
        
        figure(fH1); hold all;
        [~,ch] = megPlotMap(dataToPlot, colormapLimits, [], cmapData, {'Synch' labels{el}}, [],[], 'isolines', contourLines);
        h = findobj(gca,'Type','contour');
        h.LineWidth = lw;
        set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);
        
        % Plot contour
        figure(fH2); hold all;
        contourf(h.XData, h.YData, h.ZData, contourLines, 'LineColor',colorContourLines(el,:), 'Fill','off','LineWidth',2);
        colorbar off;
        
        %% ASYNCHRONOUS (KAPPA=0)
        figure(fH1);
        subplot(nrows, ncols, ncols+el);
        dataToPlot = squeeze(w.V1i(s,el,:))./maxSynFullStim;
        
        % Check at what point to draw contour lines
        if ~isempty(contourPercentile)
            if contourPercentile > 10 % Plot at given data percentile
                contourLines = [1 1]*prctile(dataToPlot, contourPercentile);
            else % Plot nr of contours at equo-distance percentiles
                contourLines = contourPercentile;
            end
        else, contourLines = [];
        end
        
        figure(fH1);
        [~,ch] = megPlotMap(dataToPlot, colormapLimits, [], cmapData, {'Asynch' labels{el}}, [],[], 'isolines', contourLines);
        h = findobj(gca,'Type','contour');
        h.LineWidth = lw;
        set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);
        
        % Plot contour
        figure(fH3); hold all;
        contourf(h.XData, h.YData, h.ZData, contourLines, 'LineColor',colorContourLines(el,:), 'Fill','off','LineWidth',2);
        colorbar off;
        
    end

    if saveFig
        figureDir       = fullfile(fmsRootPath,'figures', subject{s}); % Where to save images?
        if ~exist(figureDir,'dir'); mkdir(figureDir); end
        
        set(0, 'currentfigure', fH1); % save as png for memory purposes
        figurewrite(fullfile(figureDir, sprintf('SupplFigure6_varyEccenPredictions_%s_%d_%s_highResFlag%d_S%d_cumulative_normalized', area, contourPercentile, headmodelType, highResSurf, s)),0,[1 300],'.',1);
        set(0, 'currentfigure', fH2);
        figurewrite(fullfile(figureDir, sprintf('SupplFigure6_Contour_Synchr_varyEccenPredictions_%s_%d_%s_highResFlag%d_cumulative_normalized', area,contourPercentile, headmodelType, highResSurf)),[],0,'.',1);
        set(0, 'currentfigure', fH3);
        figurewrite(fullfile(figureDir, sprintf('SupplFigure6_Contour_Asynchr_varyEccenPredictions_%s_%d_%s_highResFlag%d_cumulative_normalized', area,contourPercentile, headmodelType, highResSurf)),[],0,'.',1);
    end
end



%% Take mean across subjects and plot if requested
if plotMeanSubject
    
    w.V1c_mn = squeeze(nanmean(w.V1c,1));
    w.V1i_mn = squeeze(nanmean(w.V1i,1));
    
    fH4 = figure(4); clf; set(gcf, 'Color', 'w', 'Position', [1, 1, 1680, 999]); hold all;
    fH5 = figure(5); clf; set(gcf, 'Color', 'w'); hold all; megPlotMap(zeros(1,157));  drawnow; title('Synchronous sources')
    fH6 = figure(6); clf; set(gcf, 'Color', 'w'); hold all; megPlotMap(zeros(1,157));  drawnow; title('Asynchronous sources')
    
    colormapLimits = [0 1];
    
    for el = 1:nrEccenLimits
        %% Visualize predictions SYNCHRONOUS (KAPPA=100*pi)
        
        figure(fH4);
        subplot(nrows, ncols, el);
        dataToPlot = w.V1c_mn(el,:)./maxSynFullStim;
        
        % Check at what point to draw contour lines
        if ~isempty(contourPercentile)
            if contourPercentile > 10 % Plot at given data percentile
                contourLines = [1 1]*prctile(dataToPlot, contourPercentile);
            else % Plot nr of contours at equo-distance percentiles
                contourLines = contourPercentile;
            end
        else, contourLines = []; end
        
        figure(fH4);
        [~,ch] = megPlotMap(dataToPlot, colormapLimits, [], cmapData, {'Synch' labels{el}}, [],[], 'isolines', contourLines);
        h = findobj(gca,'Type','contour');
        h.LineWidth = lw;
        set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);

        figure(fH5); hold all;
        contourf(h.XData, h.YData, h.ZData, contourLines, 'LineColor',colorContourLines(el,:), 'Fill','off','LineWidth',2);
        colorbar off;
        
        %% Visualize predictions ASYNCHRONOUS (KAPPA=0)
        
        figure(fH4); hold all;
        subplot(nrows, ncols,nrEccenLimits+el);
        dataToPlot = w.V1i_mn(el,:)./maxSynFullStim;
        
        % Check at what point to draw contour lines
        if ~isempty(contourPercentile)
            if contourPercentile > 10 % Plot at given data percentile
                contourLines = [1 1]*prctile(dataToPlot, contourPercentile);
            else % Plot nr of contours at equo-distance percentiles
                contourLines = contourPercentile;
            end
        else, contourLines = []; end
        
        figure(fH4);
        [~,ch] = megPlotMap(dataToPlot, colormapLimits, [], cmapData, {'Asynch' labels{el}}, [],[], 'isolines', contourLines);
        h = findobj(gca,'Type','contour');
        h.LineWidth = lw;
        set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);
        
        figure(fH6); hold all;
        contourf(h.XData, h.YData, h.ZData, contourLines, 'LineColor',colorContourLines(el,:), 'Fill','off','LineWidth',2);
        colorbar off;

    end
    
    if saveFig
        figureDir       = fullfile(fmsRootPath,'figures', 'average'); % Where to save images?
        if ~exist(figureDir,'dir'); mkdir(figureDir); end
        
        % Plot data and save data
        set(0, 'currentfigure', fH4); % save as png for memory purposes
        figurewrite(fullfile(figureDir, sprintf('SupplFigure6_varyEccenPredictions_%s_%d_%s_highResFlag%d_AVERAGE_cumulative_normalized', area,contourPercentile, headmodelType, highResSurf)),0,[1 300],'.',1);
        set(0, 'currentfigure', fH5);
        figurewrite(fullfile(figureDir, sprintf('SupplFigure6_Contour_Synchr_varyEccenPredictions_%s_%d_%s_highResFlag%d_AVERAGE_cumulative_normalized', area,contourPercentile, headmodelType, highResSurf)),[],0,'.',1);
        set(0, 'currentfigure', fH6);
        figurewrite(fullfile(figureDir, sprintf('SupplFigure6_Contour_Asynchr_varyEccenPredictions_%s_%d_%s_highResFlag%d_AVERAGE_cumulative_normalized', area,contourPercentile, headmodelType, highResSurf)),[],0,'.',1);
    end
    
    
end
