function makeFigure4(varargin)
%
% This is a function to make Figure 4 from the manuscript about forward
% modeling coherent and incoherent neural sources to MEG responses.
%
% This figure shows the empirical finding of two different spatial patterns
% for a stimulus-locked response and an asynchronous broadband response to
% a large field flickering (12 Hz) dartboard pattern.
%
% To runs this script, you need:
% (1) the data from the denoiseproject in the data folder of its FMS
%     code repository
%
% (2) MEG_utils toolbox added to the paths. For example:
%     tbUse('ForwardModelSynchrony');
%        or to only add the MEG_utils toolbox:
%     addpath(genpath('~/matlab/git/toolboxes/meg_utils'))
%
% INPUTS:
%   [subjectsToPlot]  :  (int) subject nr you would like to plot, default is 12
%   [plotMeanSubject] :  (bool) plot average across all 12 subjets or not?
%   [saveFig]         :  (bool) save figures or not?
%
% Example 1:
%  makeFigure4('subjectsToPlot', 1, 'plotMeanSubject', false, 'saveFig', true)
% Example 2:
%  makeFigure4('subjectsToPlot', 12, 'plotMeanSubject', false, 'saveFig', true)
% Example 3:
%  makeFigure4('subjectsToPlot', 1:12, 'plotMeanSubject', true, 'saveFig', true)

%% 0. Set up paths and define parameters
p = inputParser;
p.KeepUnmatched = true;
p.addParameter('subjectsToPlot', 12);
p.addParameter('plotMeanSubject', true, @islogical)     % Plot average across subject
p.addParameter('saveFig', true, @islogical);            % Save figures
p.addParameter('useSLIncohSpectrum', true, @islogical); % Plot SL amplitudes from incoherent spectrum (default: true)
p.addParameter('doSOIcomparison', false, @islogical);   % Compare the signal for the two types of SOI (sensors of interest, requires makeFigure1 to be executed)
p.addParameter('contourPercentile', 93.6, @isnumeric);      % Percentile of the data to draw contour lines? There are two ways to define: single int above or below 10
                                                            %   default: 97.5th percentile to select Top 10 channels. Or use 90.4 for Top 15 channels
                                                            %   alternative: get contour lines at equal percentiles of data, use any integer under 10. E.g. 3 lines -> 25 50 75th prctle
p.addParameter('maxColormapPercentile', 97.5, @isnumeric);  % At what percentile of data are we truncating colormap?
p.addParameter('signedColorbar', true, @islogical);        % Plot signed colormap (true) or only positive values (false)? 
p.addParameter('snrThresh',1, @isnumeric);                  % Threshold amplitudes by 1 SD of SNR
p.parse(varargin{:});

% Rename variables
subjectsToPlot      = p.Results.subjectsToPlot;
plotMeanSubject     = p.Results.plotMeanSubject;
saveFig             = p.Results.saveFig;
useSLIncohSpectrum  = p.Results.useSLIncohSpectrum;
doSOIcomparison     = p.Results.doSOIcomparison;
contourPercentile   = p.Results.contourPercentile;
maxColormapPercentile = p.Results.maxColormapPercentile;
signedColorbar      = p.Results.signedColorbar;
snrThresh           = p.Results.snrThresh;

% Define subjects
subject         = {'wlsubj002', ... % S1 - Full, Left, Right stim experiment
                   'wlsubj004', ... % S2 - Full, Left, Right stim experiment
                   'wlsubj005', ... % S3 - Full, Left, Right stim experiment
                   'wlsubj006', ... % S4 - Full, Left, Right stim experiment
                   'wlsubj010', ... % S5 - Full, Left, Right stim experiment
                   'wlsubj011', ... % S6 - Full, Left, Right stim experiment
                   'wlsubj048', ... % S7 - Full  stim only experiment
                   'wlsubj046', ... % S8 - Full  stim only experiment
                   'wlsubj039', ... % S9 - Full  stim only experiment
                   'wlsubj059', ... % S10 - Full  stim only experiment
                   'wlsubj067', ... % S11 - Full  stim only experiment
                   'wlsubj070'};    % S12 - Full  stim only experiment

% Set up paths
figureDir        = fullfile(fmsRootPath, 'figures'); % Where to save images?
dataDir          = fullfile(fmsRootPath, 'data');    % Where to get data?

% Number of bootstraps when comparing sensors of interest (those that fall
% within contour lines) between stimulus-locked and broadband responses
nboot = 1000;

% Preallocate space for matrices
diffFullBlankSL = NaN(length(subject),157);
diffFullBlankBB = diffFullBlankSL;

% Load all subjects when plotting the mean
if plotMeanSubject
    subjectsToLoad = 1:length(subject);
else
    subjectsToLoad = subjectsToPlot;
end

%% 1. Load subject's data

for s = subjectsToLoad
    
    % Go from subject to session nr
    switch subject{s}
        case 'wlsubj002' % S1 - Full, Left, Right stim experiment
            whichSession = 2;
        case 'wlsubj004' % S2 - Full, Left, Right stim experiment
            whichSession = 7;
        case 'wlsubj005' % S3 - Full, Left, Right stim experiment
            whichSession = 8;
        case 'wlsubj006' % S4 - Full, Left, Right stim experiment
            whichSession = 1;
        case 'wlsubj010' % S5 - Full, Left, Right stim experiment
            whichSession = 6;
        case 'wlsubj011' % S6 - Full, Left, Right stim experiment
            whichSession = 5;
        case 'wlsubj048' % S7 - Full  stim only experiment
            whichSession = 9;
        case 'wlsubj046' % S8 - Full  stim only experiment
            whichSession = 10;
        case 'wlsubj039' % S9 - Full  stim only experiment
            whichSession = 11;
        case 'wlsubj059' % S10 - Full  stim only experiment
            whichSession = 12;
        case 'wlsubj067' % S11 - Full  stim only experiment
            whichSession = 13;
        case 'wlsubj070' % S12 - Full  stim only experiment
            whichSession = 14;
    end
    
    % Get SNR data
    data = loadData(fullfile(dataDir, subject{s}),whichSession,'SNR');
    snr(s,1,:) = data{1};
    snr(s,2,:) = data{2};
    
    % Get amplitude data
    data = loadData(fullfile(dataDir, subject{s}), whichSession,'amplitudes');
    
    % Update SL amplitudes for each subject, either with coherent or incoherent spectrum
    if useSLIncohSpectrum
        ampl{s}.sl.full  = data.sl.full;
        ampl{s}.sl.blank = data.sl.blank;
        if strcmp(subject{s},'wlsubj059')
            ampl{s}.sl.full  = data.sl.full_coherent;
            ampl{s}.sl.blank = data.sl.blank_coherent;
        end
    else
        ampl{s}.sl.full  = data.sl.full_coherent;
        ampl{s}.sl.blank = data.sl.blank_coherent;
    end
    
    % Update broadband power for each subject, always from incoherent spectrum
    ampl{s}.bb.full  = data.bb.full;
    ampl{s}.bb.blank = data.bb.blank;
    
    clear data bb sl snr_sl snr_bb data;
    
    if doSOIcomparison
        % If requested, one can calculate the amplitude for the sensors
        % predicted by the model to have the highest response. We call these
        % sensors the SOI (sensors of interest) and sensor indices are saved by
        % the figurescript makeFigure4. We can compare the SOI from simulations
        % with incoherent/random or coherent/uniform time series, or their union.
        
        %  Get sensors of interest based on forward model predictions
        %  (sensorsOfInterest is boolean of 2x157)
        tmp = load(fullfile(dataDir, subject{subjectsToPlot}, sprintf('%s_sensorsOfInterestFromPrediction', subject{subjectsToPlot})));
        
        % Row 1 SOI based on single subject uniform forward model prediction
        % Row 2 SOI based on single subject random forward model prediction
        soi = logical(tmp.sensorsOfInterest); clear sensorsOfInterest; clear tmp;
        
        % Find the sensors for the single subject for uniform and random phase
        % forward model predictions
        uniformSensors = find(soi(1,:));
        randomSensors  = find(soi(2,:));
        
        % Get union and separate sensors for only uniform or random phase predictions
        singleSubject.unionUR       = intersect(uniformSensors, randomSensors);
        singleSubject.onlyUniform   = uniformSensors(~ismember(uniformSensors,singleSubject.unionUR));
        singleSubject.onlyRandom    = randomSensors(~ismember(randomSensors,singleSubject.unionUR));
        
        clear uniformSensors randomSensors;
    end



    %% 2. Get contrast between full and blank for SL and BB data

    % Take difference between mean of full and blank epochs for each subject
    % and dataset (sl or bb)    
    diffFullBlankSL(s,:) = nanmean(ampl{s}.sl.full,1) - nanmean(ampl{s}.sl.blank,1);
    diffFullBlankBB(s,:) = nanmean(ampl{s}.bb.full,1) - nanmean(ampl{s}.bb.blank,1);
    
    % there is a difference between the datasets in terms of scaling units
    % session 1-6 are in fempto Tesla  whereas 7-12 are in Tesla
    % TODO: handle this more gracefully? Maybe just mark the sessions?
    if max(diffFullBlankSL(s,:)) < 1^-14
        diffFullBlankSL(s,:) = diffFullBlankSL(s,:) .* 10^15;
        diffFullBlankBB(s,:) = diffFullBlankBB(s,:) .* 10^15 .* 10^15;
    end
end

for s = subjectsToPlot
    %% 3. Plot subject
    snrThreshMask.sl.single = abs(squeeze(snr(s,1,:))) > snrThresh;
    snrThreshMask.bb.single = abs(squeeze(snr(s,2,:))) > snrThresh;

    dataToPlot   = cat(1, diffFullBlankSL(s,:) .* snrThreshMask.sl.single', ...
        diffFullBlankBB(s,:) .* snrThreshMask.bb.single');
 
    if intersect(s, [1:7,10,12]) % [HACK: for some reason the preprocessing steps didn't take out bad channel 98 in most sessions]
         dataToPlot(:,98) = NaN; % Check if this is always a bad channel, and should be dropped earlier in the analysis
    end
     
    fig_ttl       = {sprintf('Figure4_Observed_MEG_Data_incohSpectrum%d_prctile%2.1f_S%d', useSLIncohSpectrum, contourPercentile, s), ...
                    sprintf('Figure4_Contour_incohSpectrum%d_prctile%2.1f_S%d', useSLIncohSpectrum, contourPercentile, s)};
    sub_ttl       = {sprintf('Stimulus locked S%d', s), ...
                     sprintf('Broadband S%d', s)};
    markerType    = '.';
    colorContours = {'y','b'};
    
    % Plot it!
    if saveFig
        figureDirSubj = fullfile(figureDir, subject{s});
        if ~exist(figureDirSubj, 'dir'); mkdir(figureDirSubj); end
    end
    
    visualizeSensormaps(dataToPlot, maxColormapPercentile, contourPercentile, signedColorbar, colorContours, markerType, fig_ttl, sub_ttl, saveFig, figureDirSubj);
end

%% 4. Plot average across subjects if requested

if plotMeanSubject
    
    % Get snr mask for group data
    snrThreshMask.sl.group  = abs(squeeze(nanmean(snr(:,1,:),1))) > snrThresh;
    snrThreshMask.bb.group  = abs(squeeze(nanmean(snr(:,2,:),1))) > snrThresh;
    
    % Concatenate data
    dataToPlot      = cat(1, nanmean(diffFullBlankSL,1) .* snrThreshMask.sl.group', ...
                             nanmean(diffFullBlankBB,1) .* snrThreshMask.bb.group');
    dataToPlot(:,98) = NaN;
    
    % Define figure and subfigure titles
    fig_ttl         = {sprintf('Figure4_Observed_MEG_Data_incohSpectrum%d_prctile%2.1f_AVERAGE', useSLIncohSpectrum, contourPercentile), ...
                       sprintf('Figure4_Contour_incohSpectrum%d_prctile%2.1f_AVERAGE', useSLIncohSpectrum, contourPercentile)};
    sub_ttl         = {sprintf('Stimulus locked Average N = %d', length(subject)), ...
                       sprintf('Broadband Average N = %d', length(subject))};
    
    % Make figure dir for average subject, if non-existing
    figureDirAvg       = fullfile(figureDir,'average'); % Where to save images?
    if ~exist(figureDirAvg,'dir'); mkdir(figureDirAvg); end
    
    % Plot it!
    visualizeSensormaps(dataToPlot, maxColormapPercentile, contourPercentile, signedColorbar, colorContours, markerType, fig_ttl, sub_ttl, saveFig, figureDirAvg);
end

%% 5. If requested: Make barplot of bootstrapped data of sensors that fall within contour 
% lines and compare across subjects anatomy and forward model phase

if doSOIcomparison
    
    % Bootstrap across epochs if requested to compare SOIs
    fns = fieldnames(singleSubject);
    
    % Loop over fieldnames (sensors of interest)
    for ii = 1:numel(fns)
        
        % Example subject
        bootstat{1}.sl.full.(fns{ii})  = bootstrp(nboot, @(x) nanmean(x,1), ampl{1}.sl.full(:,singleSubject.(fns{ii})));
        bootstat{1}.sl.blank.(fns{ii})  = bootstrp(nboot, @(x) nanmean(x,1), ampl{1}.sl.blank(:,singleSubject.(fns{ii})));
        
        % take difference between full and blank
        bootstat{1}.sl.diff.(fns{ii}) = mean([bootstat{1}.sl.full.(fns{ii}) - bootstat{1}.sl.blank.(fns{ii})],2);
        
        bootstat{1}.bb.full.(fns{ii})  = bootstrp(nboot, @(x) nanmean(x,1), ampl{1}.bb.full(:,singleSubject.(fns{ii})));
        bootstat{1}.bb.blank.(fns{ii})  = bootstrp(nboot, @(x) nanmean(x,1), ampl{1}.bb.blank(:,singleSubject.(fns{ii})));
        
        bootstat{1}.bb.diff.(fns{ii}) = mean([bootstat{1}.bb.full.(fns{ii}) - bootstat{1}.bb.blank.(fns{ii})],2);
        
        if plotMeanSubject
            % Get sensors of interest from forward model
            tmp = load(fullfile(dataDir, 'average', sprintf('Average_sensorsOfInterestFromPrediction_%s', area)));
            soi = logical(tmp.sensorsOfInterest); clear sensorsOfInterest; clear tmp;
            
            % And again based on the average across subjects
            uniformSensors = find(soi(1,:));
            randomSensors = find(soi(2,:));
            
            % Get union and separate sensors for only uniform or random phase predictions
            averageSubject.unionUR      = intersect(uniformSensors, randomSensors);
            averageSubject.onlyUniform  = uniformSensors(~ismember(uniformSensors,averageSubject.unionUR));
            averageSubject.onlyRandom   = randomSensors(~ismember(randomSensors,averageSubject.unionUR));
            
            % Group average
            bootstatGroup.sl.(fns{ii}) = bootstrp(nboot, @(x) nanmean(x,2), diffFullBlankSL(:,averageSubject.(fns{ii})));
            bootstatGroup.bb.(fns{ii}) = bootstrp(nboot, @(x) nanmean(x,2), diffFullBlankBB(:,averageSubject.(fns{ii})));
        end
    end
    
    % Make figure
    figure; set(gcf, 'Position', [509, 238, 672, 1072], 'Color','w')
    labels = {'Union', 'Only Uniform', 'Only Random'};
    ttls   = {'SL one subject',  'SL group average', 'BB one subject', 'BB group average'};
    
    dataAllBox = {bootstat{1}.sl.diff, bootstatGroup.sl, bootstat{1}.bb.diff, bootstatGroup.bb};
    
    ylims = [-10 50; -10 50; -2 3; -2 3];
    
    for dd = 1:4
        subplot(2,2,dd)
        
        boxplot([nanmean(dataAllBox{dd}.unionUR,2), nanmean(dataAllBox{dd}.onlyUniform,2), nanmean(dataAllBox{dd}.onlyRandom,2)], ...
            'PlotStyle','traditional', 'Widths',0.2,'MedianStyle','line','Colors','mrb'); hold on
        plot([-0.5 4.5], [0 0],'k', 'LineWidth',2)
        
        box off; ylabel('Amplitudes (fT)');
        set(gca,'TickDir','out', 'XTickLabel', labels, 'XTickLabelRotation',45, 'FontSize',12,'LineWidth',2, 'YLim', ylims(dd,:));
        title(ttls{dd});
        
    end
    
    if saveFig
        hgexport(gcf,fullfile(figureDir,'Figure4_SOI_boxplot'))
    end
    
end

return
%% OBSOLETE code

% %% Get contrast by bootstrapping over epochs
%     bootstat = cell(1,length(subject));
%
%     % Loop over fieldnames (sensors of interest)
%     for s = 1:length(subject)
%
%         % Example subject
%         bootstat{s}.sl.full.all  = bootstrp(nboot, @(x) nanmean(x,1), ampl{s}.sl.full);
%         bootstat{s}.sl.blank.all = bootstrp(nboot, @(x) nanmean(x,1), ampl{s}.sl.blank);
%
%         % take difference between full and blank
%         bootstat{s}.sl.diff.all  = mean([bootstat{s}.sl.full.all - bootstat{s}.sl.blank.all],1);
%
%         bootstat{s}.bb.full.all  = bootstrp(nboot, @(x) nanmean(x,1), ampl{s}.bb.full);
%         bootstat{s}.bb.blank.all = bootstrp(nboot, @(x) nanmean(x,1), ampl{s}.bb.blank);
%
%         bootstat{s}.bb.diff.all  = mean([bootstat{s}.bb.full.all - bootstat{s}.bb.blank.all],1);
%
%     end
%
%     if plotMeanSubject
%         % Group average
%         bootstatGroup.sl.all = bootstrp(nboot, @(x) nanmean(x,1), diffFullBlankSL);
%         bootstatGroup.bb.all = bootstrp(nboot, @(x) nanmean(x,1), diffFullBlankBB);
%     end


% Plot mean of 6 subjects

%
% % Define color range
% cLims = [-1 1]*prctile(sl_snr, 95);
%
% % Plot average subject
% fH = figure('position',[1,600,1400,800]); set(gcf, 'Name', 'Figure 2B, Average across subject', 'NumberTitle', 'off');
% [fH,ch] = megPlotMap(sl_snr,cLims,fH,'bipolar',[],data_hdr,cfg, ...
%     'isolines', contourLim*max(cLims)*[1 1], ...
%     'highlightchannel', sl_snr > contourLim*max(cLims), ...
%     'highlightmarker', '*');
%
% % colormap(bipolar);
% c = findobj(fH,'Type','Contour'); c.LineWidth = 4;
%
%
%
% set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);
%
% cLims = [-1 1]*prctile(bb_snr, 95);
%
% subplot(1,2,2)
% [~,ch] = megPlotMap(bb_snr,cLims,gcf,'bipolar',[],data_hdr,cfg,'isolines', contourLim*max(cLims)*[1 1]); colormap(bipolar);
% set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);
%
% if saveFig
%     figurewrite(fullfile(figureDir,'Figure2_SLBB_average'),[],0,'.',1);
%     hgexport(gcf,fullfile(figureDir,'Figure2_SLBB_average_2'))
% end

% (more obsolete?) Plot mean of 6 subjects with CI from contour lines

% allData  = {squeeze(sl_all)',squeeze(bb_all)'};
% clims    = {climsSL,climsBB};
% nBoot    = 50;
%
%
% meshplotFigHandleFun = @(x) megPlotMap(x, [],[],bipolar,[],[],[],'isolines', 1);
% getContourStruct     = @(x) findobj(x.Children,'Type','Contour');
% getXYZData           = @(x) cat(3,x.XData,x.YData,x.ZData);
%
% ft_warning off
%
% for d = 1:2
%
%     bootData = allData{d};
%     bootstat = bootstrp(nBoot, @(x) mean(x,1), bootData);
%
%     for ii = 1:nBoot
%         close all;
%         c              = getContourStruct(meshplotFigHandleFun(bootstat(ii,:)));
%         xyz(ii,:,:,:)= getXYZData(c);
%         clear c
%     end
%
%
%     figure('position',[1,600,1400,800]); set(gcf, 'Name', 'Figure 2B, Average across subject', 'NumberTitle', 'off');
%     subplot(121); megPlotMap(mean(bootData,1),clims{d},[],bipolar); hold all;
%     for ii = 1:nBoot
%         contour(squeeze(xyz(ii,:,:,1)),squeeze(xyz(ii,:,:,2)),squeeze(xyz(ii,:,:,3)),1, 'k-');
%     end
% %     set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);
%
%     figurewrite(fullfile(figureDir, sprintf('Figure2_dataCI_average%d',d)), [],0,'.',1);
%
% end
%
% ft_warning on


