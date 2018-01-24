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

% Set up paths
figureDir       = fullfile(fmsRootPath, 'figures'); % Where to save images?
dataDir         = fullfile(fmsRootPath, 'data');    % Where to get data?
saveFigures     = true;     % Save figures in the figure folder?
cfg             = [];
data_hdr        = [];

% Use first stimulus (full field), define how to combine data across subjects
contrasts     = [1 0 0];
contrasts     = bsxfun(@rdivide, contrasts, sqrt(sum(contrasts.^2,2)));
computeSNR    = @(x) nanmean(x,3) ./ nanstd(x, [], 3);

% What plotting range for individual example and average 
climsSLone = [-1 1] * 20;
climsBBone = [-1 1] * 5;

climsSLave = [-1 1] * 15;
climsBBave = [-1 1] * 2;

% Predefine tickmark position for colorbar
yscaleAB = [-6,-3,0,3,6];


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


%% 2. Plot one subject

% Plot subject 1 (wl_subj002)
figure('position',[1,600,1400,800]); set(gcf, 'Name', 'Figure 2A, Example subject', 'NumberTitle', 'off');
subplot(1,2,1)
[~,ch] = megPlotMap(sl_all(:,:,1),climsSLone,gcf,'bipolar',[],data_hdr,cfg,'isolines', 0.64*max(climsSLone)*[1 1]); colormap(bipolar);
set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);

subplot(1,2,2)
[~,ch] = megPlotMap(bb_all(:,:,1),climsBBone,gcf,'bipolar',[],data_hdr,cfg,'isolines', 0.32*max(climsBBone)*[1 1]); colormap(bipolar);
set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);

if saveFigures
     figurewrite(fullfile(figureDir,'Figure2_SLBB_onesubject'),[],0,'.',1);
        hgexport(gcf,fullfile(figureDir,'Figure2_SLBB_onesubject_2'))

end

%% 3. Plot mean of 6 subjects

% Get stimulus-locked snr across subjects
sl_snr = nanmean(sl_all,3);
% Get broadband snr (after denoising) across subjects
bb_snr = nanmean(bb_all,3);

% Define nans
idx = isfinite(sl_snr);

% Plot average subject
figure('position',[1,600,1400,800]); set(gcf, 'Name', 'Figure 2B, Average across subject', 'NumberTitle', 'off');
subplot(1,2,1)
[~,ch] = megPlotMap(sl_snr,climsSLave,gcf,'bipolar',[],data_hdr,cfg, 'isolines', 0.64*max(climsSLave)*[1 1]); colormap(bipolar);
set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);
subplot(1,2,2)
[~,ch] = megPlotMap(bb_snr,climsBBave,gcf,'bipolar',[],data_hdr,cfg,'isolines', 0.30*max(climsBBave)*[1 1]); colormap(bipolar);
set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);x

if saveFigures
     figurewrite(fullfile(figureDir,'Figure2_SLBB_average'),[],0,'.',1);
    hgexport(gcf,fullfile(figureDir,'Figure2_SLBB_average_2'))
end

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

return
