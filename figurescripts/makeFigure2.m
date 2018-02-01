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
subject = {'wl_subj002','wl_subj004','wl_subj005','wl_subj006','wl_subj010','wl_subj011'};

% Which example subject to show?
exampleSubject = 1;

% What's the plotting range for individual example and average across
% subjects?
contourmapPercentile   = 93.6; % draw contour line at what fraction of the colormap?
colormapPercentile     = 97.5; % percentile of data to use for max/min limits of colorbar

% middleSensors = [13 23, 2, 20, 25, 43];

% Set up paths
figureDir       = fullfile(fmsRootPath, 'figures'); % Where to save images?
dataDir         = fullfile(fmsRootPath, 'data');    % Where to get data?
saveFigures     = true;     % Save figures in the figure folder?

% Use first stimulus (full field), define how to combine data across subjects
contrasts     = [1 0 0];
contrasts     = bsxfun(@rdivide, contrasts, sqrt(sum(contrasts.^2,2)));
computeSNR    = @(x) nanmean(x,3) ./ nanstd(x, [], 3);

% Predefine tickmark position for colorbar
% yscaleAB = [-6,-3,0,3,6];


%% 1. Load subject's data

for s = 1:length(subject)
    
    switch subject{s}
        % Go from subject to session nr
        case 'wl_subj002'
            whichSubject = 2;
        case 'wl_subj004'
            whichSubject = 7;
        case 'wl_subj005'
            whichSubject = 8;
        case 'wl_subj006'
            whichSubject = 1;
        case 'wl_subj010'
            whichSubject = 6;
        case 'wl_subj011'
            whichSubject = 5;
    end
    
    data = loadData(dataDir,whichSubject);
    bb = data{1}; sl = data{2}; clear data;
    
    % Get information, assuming this is the same for SL and BB data
    numChannels  = size(sl.results.origmodel.beta,2);
    numBoots     = size(sl.results.origmodel.beta,3);
    numContrasts = 1;
    
    % Stimulus-locked: Compute SNR for contrasts
    tmp_data = reshape(sl.results.origmodel.beta,3,[]);
    tmp      = contrasts*tmp_data;
    tmp      = reshape(tmp, numContrasts, numChannels, numBoots);
    sSL      = computeSNR(tmp)';
    
    % Broadband before denoising: Compute SNR for contrasts
    tmp_data = reshape(bb.results.origmodel.beta,3,[]);
    tmp = contrasts*tmp_data;
    tmp = reshape(tmp, numContrasts, numChannels,numBoots);
    sBB = computeSNR(tmp)';
    
    % Prepare array
    if s == 1
        sl_all = NaN(size(contrasts,1),length(sl(1).badChannels), length(subject));
        bb_all = sl_all;
    end
    
    % Update array with data converted to channel space
    sl_all(:,:,s) = to157chan(sSL', ~sl.badChannels,'nans');
    bb_all(:,:,s) = to157chan(sBB', ~bb.badChannels,'nans');
    
end

%% Get stimulus-locked snr across subjects
sl_snr = nanmean(sl_all,3);
% Get broadband snr (after denoising) across subjects
bb_snr = nanmean(bb_all,3);


%% 2. Plot one subject and average across subjects

dataAll      = cat(1,sl_all(:,:,exampleSubject), bb_all(:,:,exampleSubject), sl_snr, bb_snr);
colorMarkers = {'r','b', 'r', 'b'};
fig_ttl      = {'Figure2_Data', 'Figure2_Sl_and_Broadband_Compared'};
sub_ttl          = {'Stimulus Locked S1', ...
                'Broadband S2', ...
                'Stimulus Locked S1-S6', ...
                'Broadband S1-S6'};
markerType   = '.';

visualizeSensormaps(dataAll, colormapPercentile, contourmapPercentile, colorMarkers, markerType, fig_ttl, sub_ttl, saveFigures, figureDir)


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


