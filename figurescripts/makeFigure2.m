function makeFigure2(exampleSubject)

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

% Which subjects to average?
%   Full  only: 'wlsubj048', 'wlsubj046','wl_subj039','wl_subj059', 'wl_subj067'
%   Full, Left, Right: 'wl_subj002','wl_subj004','wl_subj005','wl_subj006','wl_subj010','wl_subj011'
subject         = {'wlsubj070'};%{'wl_subj002','wl_subj004','wl_subj005','wl_subj006','wl_subj010','wl_subj011','wlsubj048', 'wlsubj046','wl_subj039','wl_subj059', 'wl_subj067'};
if nargin < 1; exampleSubject  = 1; end % Which example subject to show if not defined

% Set up paths
figureDir              = fullfile(fmsRootPath, 'figures', subject{exampleSubject}); % Where to save images?
dataDir                = fullfile(fmsRootPath, 'data');    % Where to get data?
saveFigures            = true;      % Save figures in the figure folder?
plotMeanSubject        = false;     % Plot average subject?
useSLIncohSpectrum     = true;      % Plot SL amplitudes from incoherent spectrum (default: true)
doSOIcomparison        = false;     % Compare the signal for the two types of SOI (sensors of interest, requires makeFigure1 to be executed)

% What contour and color map percentile to use to define the limits
contourmapPercentile   = 93.6; % draw contour line at what fraction of the colormap?
colormapPercentile     = 97.5; % percentile of data to use for max/min limits of colorbar
snrThresh              = 1;    % Threshold amplitudes by 1 SD of SNR

% Number of bootstraps
nboot = 1000;

% Predefine tickmark position for colorbar
% yscaleAB = [-6,-3,0,3,6];


%% 1. Load subject's data

for s = 1:length(subject)
    switch subject{s}
        % Go from subject to session nr
        case 'wl_subj002'
            whichSession = 2;
        case 'wl_subj004'
            whichSession = 7;
        case 'wl_subj005'
            whichSession = 8;
        case 'wl_subj006'
            whichSession = 1;
        case 'wl_subj010'
            whichSession = 6;
        case 'wl_subj011'
            whichSession = 5;
        case 'wlsubj048'
            whichSession = 9; % Full field Only
        case 'wlsubj046'
            whichSession = 10; % Full field Only
        case 'wl_subj039'
            whichSession = 11; % Full field Only
        case 'wl_subj059'
            whichSession = 12; % Full field Only
        case 'wl_subj067'
            whichSession = 13; % Full field Only
        case 'wlsubj070'
            whichSession = 14; % Full field Only
    end
    
    
    
    % Get SNR data
    data = loadData(fullfile(dataDir, subject{s}),whichSession,'SNR');
    snr(s,1,:) = data{1};
    snr(s,2,:) = data{2};
   
    % Get amplitude data
    data = loadData(fullfile(dataDir, subject{s}), whichSession,'amplitudes');
    
    % Update array with data converted to channel space
    if useSLIncohSpectrum
        ampl{s}.sl.full  = data.sl.full;
        ampl{s}.sl.blank = data.sl.blank;
        
    else
        ampl{s}.sl.full  = data.sl.full_coherent;
        ampl{s}.sl.blank = data.sl.blank_coherent;       
    end
    

    ampl{s}.bb.full  = data.bb.full;
    ampl{s}.bb.blank = data.bb.blank;

    clear data bb sl snr_sl snr_bb data;
    
    if doSOIcomparison
        
        %  Get sensors of interest based on forward model predictions
        %  (sensorsOfInterest is booloean of 2x157)

        % Row 1 SOI based on single subject uniform forward model prediction,
        % Row 2 SOI based on single subject random forward model prediction,
        tmp = load(fullfile(dataDir, subject{exampleSubject}, sprintf('%s_sensorsOfInterestFromPrediction', subject{exampleSubject})));
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
end


%% Get contrast between full and blank for SL and BB data

% Init matrices
diffFullBlankSL = NaN(length(subject),size(ampl{s}.sl.full,2));
diffFullBlankBB = diffFullBlankSL;

% Take difference between mean of full and blank epochs for each subject
% and dataset (sl or bb)
for s = 1:length(subject)
    diffFullBlankSL(s,:) = nanmean(ampl{s}.sl.full,1) - nanmean(ampl{s}.sl.blank,1);
    diffFullBlankBB(s,:) = nanmean(ampl{s}.bb.full,1) - nanmean(ampl{s}.bb.blank,1);
end


%% Bootstrap across epochs if requested

if doSOIcomparison 

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
            % Group average
            bootstatGroup.sl.(fns{ii}) = bootstrp(nboot, @(x) nanmean(x,2), diffFullBlankSL(:,averageSubject.(fns{ii})));
            bootstatGroup.bb.(fns{ii}) = bootstrp(nboot, @(x) nanmean(x,2), diffFullBlankBB(:,averageSubject.(fns{ii})));
        end

    end
    
end


%% 2. Plot one subject

snrThreshMask.sl.single = abs(squeeze(snr(exampleSubject,1,:))) > snrThresh;
snrThreshMask.bb.single = abs(squeeze(snr(exampleSubject,2,:))) > snrThresh;

dataToPlot   = cat(1, diffFullBlankSL(1,:) .* snrThreshMask.sl.single', ...
                 diffFullBlankBB(1,:) .* snrThreshMask.bb.single');
% dataToPlot(:,98)=NaN; % Check if this is always a bad channel, or only
% for a certain subject
fig_ttl      = {sprintf('Figure2_Observed_MEG_Data_%d', useSLIncohSpectrum), sprintf('Figure2_Sl_and_Broadband_Compared_%d', useSLIncohSpectrum)};
sub_ttl      = {sprintf('Stimulus locked S%d', exampleSubject), ...
                sprintf('Broadband S%d', exampleSubject)};

% Plot it!
visualizeSensormaps(dataToPlot, colormapPercentile, contourmapPercentile, [], [], fig_ttl, sub_ttl, saveFigures, figureDir);


%% 3. Plot average across subjects if requested

if plotMeanSubject
    
    % Get sensors of interest from forward model
    tmp = load(fullfile(dataDir, subject{exampleSubject}, 'Average_sensorsOfInterestFromPrediction'));
    soi = logical(tmp.sensorsOfInterest); clear sensorsOfInterest; clear tmp;
    
    % And again based on the average across subjects
    uniformSensors = find(soi(1,:));
    randomSensors = find(soi(2,:));
    
    % Get union and separate sensors for only uniform or random phase predictions
    averageSubject.unionUR      = intersect(uniformSensors, randomSensors);
    averageSubject.onlyUniform  = uniformSensors(~ismember(uniformSensors,averageSubject.unionUR));
    averageSubject.onlyRandom   = randomSensors(~ismember(randomSensors,averageSubject.unionUR));

    % Get snr mask for group data    
    snrThreshMask.sl.group  = abs(squeeze(nanmean(snr(:,1,:),1))) > snrThresh;
    snrThreshMask.bb.group  = abs(squeeze(nanmean(snr(:,2,:),1))) > snrThresh;
    
    % Concatenate data
    dataToPlot      = cat(1, nanmean(diffFullBlankSL,1) .* snrThreshMask.sl.group', ...
                      nanmean(diffFullBlankBB,1) .* snrThreshMask.bb.group');
                  
    % Define figure and subfigure titles  
    fig_ttl         = {'Figure2_Observed_MEG_Data', 'Figure2_Sl_and_Broadband_Compared'};
    sub_ttl         = {sprintf('Stimulus locked Average N = %d', length(subject)), ...
                        sprintf('Broadband Average N = %d', length(subject))};
    
    % Make figure dir for average subject, if non-existing             
    if ~exist(fullfile(fmsRootPath, 'figures', 'average'),'dir'); mkdir(fullfile(fmsRootPath, 'figures', 'average')); end
    figureDir       = fullfile(fmsRootPath,'figures', 'average'); % Where to save images?

    % Plot it!
    visualizeSensormaps(dataToPlot, colormapPercentile, contourmapPercentile, [], [], fig_ttl, sub_ttl, saveFigures, figureDir);


    %% 4. Make barplot of bootstrapped data of sensors that fall within contour lines and compare across subjects anatomy and forward model phase
    if doSOIcomparison
        
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
end

return
%% OBSOLETE 3. Plot mean of 6 subjects
%
%
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

%% (obsolete?) Plot mean of 6 subjects with CI from contour lines

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


