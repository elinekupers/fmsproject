function makeFigure2(subjectsToPlot)

% This is a function to make Figure 2 from the manuscript about forward
% modeling coherent and incoherent neural sources to MEG responses.

% This figure shows the empirical finding of two different spatial patterns
% for a stimulus-locked response and an asynchronous broadband response to
% a large field flickering (12 Hz) dartboard pattern.

% To runs this script, you need:
% (1) the data from the denoiseproject in the data folder of its FMS
%     code repository
%
% (2) MEG_utils toolbox added to the paths. For example:
%     tbUse('ForwardModelSynchrony');
%        or to only add the MEG_utils toolbox:
%     addpath(genpath('~/matlab/git/toolboxes/meg_utils'))


%% 0. Set up paths and define parameters

% Which example subject to show if not defined
if nargin < 1; subjectsToPlot  = 12; end 

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
figureDir              = fullfile(fmsRootPath, 'figures'); % Where to save images?
dataDir                = fullfile(fmsRootPath, 'data');    % Where to get data?
saveFigures            = true;      % Save figures in the figure folder?
plotMeanSubject        = false;     % Plot average subject?
useSLIncohSpectrum     = false;      % Plot SL amplitudes from incoherent spectrum (default: true)
doSOIcomparison        = false;     % Compare the signal for the two types of SOI (sensors of interest, requires makeFigure1 to be executed)

% Two ways of plotting contours
% (1) one contour line at the percentile of data (say 90.4/100).
% (2) number of contour lines, dividing data into equal groups (use one number under 10)
%   for example, if contourPercentile=3, you draw 3 lines at the 25, 50 and 75th percentile
contourPercentile     = 93.6;
% contourPercentile      = 3;
colormapPercentile     = 97.5; % percentile of data to use for max/min limits of colorbar
snrThresh              = 0;    % Threshold amplitudes by 1 SD of SNR

% Number of bootstraps
nboot = 1000;

% Predefine tickmark position for colorbar
% yscaleAB = [-6,-3,0,3,6];

% Preallocate space for matrices
diffFullBlankSL = NaN(length(subject),157);
diffFullBlankBB = diffFullBlankSL;

%% 1. Load subject's data

for s = subjectsToPlot
    
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
        % the figurescript makeFigure1. We can compare the SOI from simulations
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

    %% 3. Plot example subject
    snrThreshMask.sl.single = abs(squeeze(snr(s,1,:))) > snrThresh;
    snrThreshMask.bb.single = abs(squeeze(snr(s,2,:))) > snrThresh;

    dataToPlot   = cat(1, diffFullBlankSL(s,:) .* snrThreshMask.sl.single', ...
        diffFullBlankBB(s,:) .* snrThreshMask.bb.single');
    dataToPlot(:,98) = NaN; % Check if this is always a bad channel, or only
    % for a certain subject
    fig_ttl      = {sprintf('Figure2_Observed_MEG_Data_incohSpectrum%d_contour%2.1f', useSLIncohSpectrum, contourPercentile), sprintf('Figure2_Sl_and_Broadband_Compared_incohSpectrum%d_contour%2.1f', useSLIncohSpectrum, contourPercentile)};
    sub_ttl      = {sprintf('Stimulus locked S%d', s), ...
        sprintf('Broadband S%d', s)};

    % Plot it!
    if saveFigures
        figureDir = fullfile(figureDir, subject{s});
        if ~exist(figureDir, 'dir'); mkdir(figureDir); end
    end
    
    visualizeSensormaps(dataToPlot, colormapPercentile, contourPercentile, [], [], fig_ttl, sub_ttl, saveFigures, figureDir);
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
    fig_ttl         = {sprintf('Figure2_Observed_MEG_Data_incohSpectrum%d_contour%2.1f', useSLIncohSpectrum, contourPercentile), sprintf('Figure2_Sl_and_Broadband_Compared_incohSpectrum%d_contour%2.1f', useSLIncohSpectrum, contourPercentile)};
    sub_ttl         = {sprintf('Stimulus locked Average N = %d', length(subject)), ...
        sprintf('Broadband Average N = %d', length(subject))};
    
    % Make figure dir for average subject, if non-existing
    if ~exist(fullfile(fmsRootPath, 'figures', 'average'),'dir'); mkdir(fullfile(fmsRootPath, 'figures', 'average')); end
    figureDir       = fullfile(fmsRootPath,'figures', 'average'); % Where to save images?
    
    % Plot it!
    visualizeSensormaps(dataToPlot, colormapPercentile, contourPercentile, [], [], fig_ttl, sub_ttl, saveFigures, figureDir);
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
    
    if saveFigures
        hgexport(gcf,fullfile(figureDir,'Figure2_SOI_boxplot'))
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
% if saveFigures
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


