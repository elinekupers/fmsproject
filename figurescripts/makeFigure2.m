function makeFigure2()

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
subject         = {'wlsubj048', 'wlsubj046','wl_subj039','wl_subj059', 'wl_subj067'};

% Which example subject to show?
exampleSubject = 1;

% Set up paths
figureDirSub    = fullfile(fmsRootPath, 'figures', subject{exampleSubject}); % Where to save images?
figureDir       = fullfile(fmsRootPath, 'figures'); % Where to save images?
dataDir         = fullfile(fmsRootPath, 'data');    % Where to get data?
saveFigures     = true;     % Save figures in the figure folder?

% What's the plotting range for individual example and average across
% subjects?
tmp = load(fullfile(dataDir, subject{exampleSubject}, sprintf('%s_prediction', subject{exampleSubject})));
contourmapPercentile = tmp.dataAll;

% contourmapPercentile   = 93.6; % draw contour line at what fraction of the colormap? 
colormapPercentile     = 97.5; % percentile of data to use for max/min limits of colorbar
snrThresh              = 1;    % Threshold amplitudes by 1 SD of SNR

% Number of bootstraps
nboot = 1000;

% Predefine tickmark position for colorbar
% yscaleAB = [-6,-3,0,3,6];


%% 1. Load subject's data

% predefine cells
soi_full_sl =  cell(length(subject),2);
soi_blank_sl = cell(length(subject),2);
soi_full_bb =  cell(length(subject),2);
soi_blank_bb = cell(length(subject),2);

% Get sensors of interest based on forward model predictions
% sensorsOfInterest is booloean of 4x157

% Row 1 SOI based on single subject uniform forward model prediction,
% Row 2 SOI based on average uniform forward model prediction,
% Row 3 SOI based on single subject random forward model prediction,
% Row 4 SOI based on average subject random forward model prediction
tmp = load(fullfile(dataDir, subject{exampleSubject}, sprintf('%s_sensorsOfInterestFromPrediction', subject{exampleSubject})));
soi = logical(tmp.sensorsOfInterest); clear sensorsOfInterest; clear tmp;

% Find the sensors for the single subject for uniform and random phase
% forward model predictions
uniformSensors = find(soi(1,:));
randomSensors = find(soi(3,:));

% Get union and separate sensors for only uniform or random phase predictions
singleSubject.unionUR       = intersect(uniformSensors, randomSensors);
singleSubject.onlyUniform   = uniformSensors(~ismember(uniformSensors,singleSubject.unionUR));
singleSubject.onlyRandom    = randomSensors(~ismember(randomSensors,singleSubject.unionUR));

% And again based on the average across subjects
uniformSensors = find(soi(2,:));
randomSensors = find(soi(4,:));

% Get union and separate sensors for only uniform or random phase predictions
averageSubject.unionUR      = intersect(uniformSensors, randomSensors);
averageSubject.onlyUniform  = uniformSensors(~ismember(uniformSensors,averageSubject.unionUR));
averageSubject.onlyRandom   = randomSensors(~ismember(randomSensors,averageSubject.unionUR));

clear uniformSensors randomSensors;



% Get data by looping over subjects
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
    end
    
    
    
    % Get SNR data
    data = loadData(fullfile(dataDir, subject{s}),whichSession,'SNR');
    bb = data{1};
    sl = data{2};

    
    % Pick full field condition (first one, or the only one)
    if whichSession >8; whichCondition = 1; else whichCondition = [1 0 0]; end
    
        
    % get stimulus-locked snr
    snr_sl = getsignalnoise(sl.results.origmodel(1), whichCondition, 'SNR',sl.badChannels);
    
    % get broadband snr for before
    snr_bb = getsignalnoise(bb.results.origmodel(1), whichCondition, 'SNR',bb.badChannels);
    
    % Account for NaNs in the data
    snr(s,1,:) = to157chan(snr_sl,~sl.badChannels,'nans');
    snr(s,2,:) = to157chan(snr_bb,~bb.badChannels,'nans');  
    
    % coherent spectrum
    data = loadData(fullfile(dataDir, subject{s}),whichSession,'amplitudes');
    bb = data.bb;
    sl = data.sl;
%     slCoh = getstimlocked_coherent(
    


    
   
    % Get amplitude data
    [data, badChannels] = loadData(fullfile(dataDir, subject{s}), whichSession,'amplitudes');
    
    % Update array with data converted to channel space
    ampl{s}.sl.full = to157chan(data.sl.full, ~badChannels,'nans');
    ampl{s}.sl.blank = to157chan(data.sl.blank, ~badChannels,'nans');
    
    ampl{s}.bb.full = to157chan(data.bb.full, ~badChannels,'nans');
    ampl{s}.bb.blank = to157chan(data.bb.blank, ~badChannels,'nans');
    
    clear bb sl snr_sl snr_bb data;
    
end


%% Bootstrap across epochs


fns = fieldnames(singleSubject);

% Init matrices
diffFullBlankSL = NaN(length(subject),size(ampl{s}.sl.full,2));
diffFullBlankBB = diffFullBlankSL;

% Take difference between mean of full and blank epochs for each subject
% and dataset (sl or bb)
for s = 1:length(subject)
    diffFullBlankSL(s,:) = nanmean(ampl{s}.sl.full,1) - nanmean(ampl{s}.sl.blank,1);
    diffFullBlankBB(s,:) = nanmean(ampl{s}.bb.full,1) - nanmean(ampl{s}.bb.blank,1);
end


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
        
    % Group average
    bootstatGroup.sl.(fns{ii}) = bootstrp(nboot, @(x) nanmean(x,2), diffFullBlankSL(:,averageSubject.(fns{ii})));
    bootstatGroup.bb.(fns{ii}) = bootstrp(nboot, @(x) nanmean(x,2), diffFullBlankBB(:,averageSubject.(fns{ii})));
    
end


%% 2. Plot one subject and average across subjects

snrThreshMask.sl.single = abs(squeeze(snr(exampleSubject,1,:))) > snrThresh;
snrThreshMask.sl.group  = abs(squeeze(nanmean(snr(:,1,:),1))) > snrThresh;

snrThreshMask.bb.single = abs(squeeze(snr(exampleSubject,2,:))) > snrThresh;
snrThreshMask.bb.group  = abs(squeeze(nanmean(snr(:,2,:),1))) > snrThresh;


dataAllMesh      = cat(1, diffFullBlankSL(1,:) .* snrThreshMask.sl.single', ...
                      diffFullBlankBB(1,:) .* snrThreshMask.bb.single', ...
                      nanmean(diffFullBlankSL,1) .* snrThreshMask.sl.group', ...
                      nanmean(diffFullBlankBB,1) .* snrThreshMask.bb.group');
                  
fig_ttl      = {'Figure2_Observed_MEG_Data', 'Figure2_Sl_and_Broadband_Compared'};
sub_ttl          = {sprintf('Stimulus locked S%d', exampleSubject), ...
                sprintf('Broadband S%d', exampleSubject), ...
                'Stimulus locked group average', ...
                'Broadband group average'};

visualizeSensormaps(dataAllMesh, colormapPercentile, contourmapPercentile, [], [], fig_ttl, sub_ttl, saveFigures, figureDirSub);

%% 3. Make barplot with sensors that fall within contour lines

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


