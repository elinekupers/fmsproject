function makeFigure2()

% addpath(genpath('~/matlab/git/denoiseproject'))
% addpath(genpath('~/matlab/git/toolboxes/meg_utils'))

% Which subjects to average?
subject = {'wl_subj002','wl_subj004','wl_subj005','wl_subj006','wl_subj010','wl_subj011'};

% Set up paths
figureDir       = fullfile(fmsRootPath, 'figures'); % Where to save images?
dataDir         = fullfile(dfdRootPath, 'analysis', 'data');    % Where to save data?
saveFigures     = true;     % Save figures in the figure folder?
threshold       = 0;        % Set threshold for colormap. If no threshold set value to 0
cfg             = [];
data_hdr        = [];

% What stimulus, how to combine data across subjects
contrastNames = {'Stim Full'}; %,'Stim Left','Stim Right','Left minus Right'
contrasts = [1 0 0]; %[eye(3); 0 1 -1];
contrasts = bsxfun(@rdivide, contrasts, sqrt(sum(contrasts.^2,2)));
computeSNR    = @(x) nanmean(x,3) ./ nanstd(x, [], 3);

for s = 1:length(subject)
    %     whichSubject = s;
    % Go from subject to session nr
    if strcmp(subject{s},'wl_subj002')
        whichSubject = 2;
    elseif strcmp(subject{s},'wl_subj004')
        whichSubject = 7;
    elseif strcmp(subject{s},'wl_subj005')
        whichSubject = 8;
    elseif strcmp(subject{s},'wl_subj006')
        whichSubject = 1;
    elseif strcmp(subject{s},'wl_subj010')
        whichSubject = 6;
    elseif strcmp(subject{s},'wl_subj011')
        whichSubject = 5;
    end
    
    data = prepareData(dataDir,whichSubject,5);
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


%% Plot one subject


yscaleAB = [repmat([-8,-4,0,4,8],3,1);[-5,-2.5,0,2.5,5]];
climsSL = [-25.6723,25.6723];
climsBB = [-8.4445, 8.4445];

% Subject 1 (wl_subj002)
figure('position',[1,600,1400,800]); set(gcf, 'Name', 'Figure 2A, Example subject', 'NumberTitle', 'off');
subplot(1,2,1)
[~,ch] = megPlotMap(sl_all(:,:,1),climsSL,gcf,'bipolar',[],data_hdr,cfg,'isolines', 1); colormap(bipolar);
set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);

subplot(1,2,2)
[~,ch] = megPlotMap(bb_all(:,:,1),climsBB,gcf,'bipolar',[],data_hdr,cfg,'isolines', 1); colormap(bipolar);
set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);

if saveFigures
    figurewrite(fullfile(figureDir,'Figure2_SLBB_average'),[],0,'.',1);
end

%% Plot mean of 6 subjects


% get stimulus-locked snr
sl_snr = nanmean(sl_all,3);
% get broadband snr for before and after denoising
bb_snr = nanmean(bb_all,3);

% Define nans
idx = isfinite(sl_snr);

figure('position',[1,600,1400,800]); set(gcf, 'Name', 'Figure 2B, Average across subject', 'NumberTitle', 'off');
subplot(1,2,1)
[~,ch] = megPlotMap(sl_snr,climsSL,gcf,'bipolar',[],data_hdr,cfg, 'isolines', 1); colormap(bipolar);
set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);
subplot(2,2,2)
[~,ch] = megPlotMap(bb_snr,climsBB,gcf,'bipolar',[],data_hdr,cfg,'isolines', 1); colormap(bipolar);
set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);

if saveFigures
    figurewrite(fullfile(figureDir,'Figure2_SLBB_average'),[],0,'.',1);
end

%% Plot mean of 6 subjects with CI from contour lines

bootData = {sl_snr,bb_snr};
clims = {climsSL,climsBB};
for d = 1:2
    
    
    bootData = bootData{d};
    nBoot    = 100;
    
    meshplotFigHandleFun = @(x) megPlotMap(x, [],[],bipolar,[],[],[],'isolines', 1);
    getContourStruct     = @(x) findobj(x.Children,'Type','Contour');
    getXYZData           = @(x) cat(3,x.XData,x.YData,x.ZData);
    
    bootstat = bootstrp(nBoot, @(x) nanmean(x,1), bootData);
    
    for ii = 1:nBoot
        close all;
        c              = getContourStruct(meshplotFigHandleFun(bootstat(ii,:)));
        xyz(ii,:,:,:)= getXYZData(c);
        clear c
    end
    
    figure('position',[1,600,1400,800]); set(gcf, 'Name', 'Figure 2B, Average across subject', 'NumberTitle', 'off');
   
    subplot(121); megPlotMap(bootData,clims{d},[],bipolar); hold all;
    for ii = 1:nBoot
        contour(squeeze(xyz(ii,:,:,1)),squeeze(xyz(ii,:,:,2)),squeeze(xyz(ii,:,:,3)),1, 'k-'); colorbar;
    end
    set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);

    figurewrite(fullfile(figureDir, sprintf('Figure2_dataCI_average%d',d)), [],0,'.',1);
    
end


