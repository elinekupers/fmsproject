function makeFigure3()

% This is a function to make Figure 3 from the manuscript about forward
% modeling coherent and incoherent neural sources to MEG responses.

% This figure shows overlap between model and data by plotting the contour 
% lines of the MEG forward model based on coherent and incoherent predictions 
% coming from vertices located in V1, on top of the stimulus-locked data 

% To runs this script, you need: 
% (1) the data from the denoiseproject in the data folder of its FMS 
%     code repository
% (2) Access to the SSMEG folder in the brainstorm data base    
% (3) MEG_utils and Fieldtrip toolbox added to the paths. For example:
%     tbUse('ForwardModelSynchrony');
%        or to only add the MEG_utils toolbox:
%     addpath(genpath('~/matlab/git/toolboxes/meg_utils'))

%% 0. Define paths and variables
bsDB            = '/Volumes/server/Projects/MEG/brainstorm_db/';
figureDir       = fullfile(fmsRootPath, 'figures'); % Where to save images?
saveFigures     = true;     % Save figures in the figure folder?
dataDir         = fullfile(fmsRootPath, 'data');    % Where to get data?

% Define project name, subject and data/anatomy folders
projectName    = 'SSMEG';

% Which subjects to average?
%   Full  only: 'wlsubj048', 'wlsubj046','wl_subj039','wl_subj059', 'wl_subj067'
%   Full, Left, Right: 'wl_subj002','wl_subj004','wl_subj005','wl_subj006','wl_subj010','wl_subj011'
subject         = {'wl_subj002','wl_subj004','wl_subj005','wl_subj006','wl_subj010','wl_subj011'};

% What's the plotting range
climsSL = [-20,20];
climsBB = [-6,6];
climsPred = [-1 1]*1E-4;

% Number of iterations for the random coherence prediction of the forward
% model
n        = 1000;     % number of timepoints (ms)
nrEpochs = 1;        % number of epochs

% Define vector that can truncate number of sensors 
keep_sensors = logical([ones(157,1); zeros(192-157,1)]); % Note: Figure out a more generic way to define keep_sensors

% Predefine figures
fH1 = figure(1); clf; set(fH1, 'Position', [1 1 1600 800], 'Name','Figure 3A, Data against model predictions V1 - matched');
fH3 = figure(3); clf; set(fH3, 'Position', [1 1 1600 800], 'Name','Figure 3B, Data against model predictions V1 - unmatched');

for s = 1:length(subject)
    
    d = dir(fullfile(bsDB, projectName, 'data', subject{s}, 'R*'));
    if strcmp(subject{s},'wl_subj059'); idx =2; else idx = 1; end
    BSdataDir = fullfile(d(idx).folder, d(idx).name);    
    BSanatDir = fullfile(bsDB, projectName, 'anat', subject{s});
    
    %% 1. Load relevant matrices
    
    G_constrained = getGainMatrix(BSdataDir, keep_sensors);

    % Get V1 template limited to 11 degrees eccentricity
    template = getTemplate(BSanatDir, 'V1', 11);

    % Simulate coherent and incoherent source time series and compute
    % predictions from forward model (w)
    tmp = getForwardModelPredictions(G_constrained, template.V1StimEccen, [], n, nrEpochs);
   
    % Take mean across epochs
    w.V1c(s,:) = mean(abs(tmp.c),2);
    w.V1i(s,:) = mean(abs(tmp.i),2);
    
    %% 3. Data
    
    % Get stimulus locked and broadband response
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
    
    % Load denoised data of example subject
    data = loadData(fullfile(dataDir, subject{s}),whichSession, 'SNR');
    bb = data{1};
    sl = data{2};
    
    % Pick full field condition (first one, or the only one)
    if whichSession >8; whichCondition = 1; else whichCondition = [1 0 0]; end
    
    % get stimulus-locked snr
    sl_signal = getsignalnoise(sl.results.origmodel(1),whichCondition, 'SNR',sl.badChannels);
    % get broadband snr for before and after denoising
    bb_signal = getsignalnoise(bb.results.origmodel(1),  whichCondition, 'SNR',bb.badChannels);
    
    % Account for NaNs in the data
    sl_signal = to157chan(sl_signal,~sl.badChannels,'nans');
    bb_signal = to157chan(bb_signal,~bb.badChannels,'nans');
    
    
    %% 4. Plotting to get contour lines
    figure(2); clf;
    ax1 =subplot(211);
    megPlotMap(abs(w.V1c(s,:)),climsPred,[],bipolar,[],[],[],'isolines', 3);
    c1 = findobj(ax1.Children,'Type','Contour');
    
%     ax2 =subplot(212);
%     megPlotMap(abs(w.V1i(s,:)),0.5*climsPred,[],bipolar,[],[],[],'isolines', 3);
%     c2 = findobj(ax2.Children,'Type','Contour');
    
    figure(1);
    subplot(2,length(subject),s)
    [~,ch] = megPlotMap(sl_signal,climsSL,gcf,'bipolar');
    hold on; contour(c1.XData,c1.YData, c1.ZData,3, 'k-'); colormap(bipolar);
    set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);
    
%     subplot(2,length(subject),s+length(subject))
%     [~,ch] = megPlotMap(bb_signal,climsBB,gcf,'bipolar');
%     hold on; contour(c2.XData,c2.YData, c2.ZData,3, 'k-'); colormap(bipolar);
%     set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);
%     
%     figure(3);
%     subplot(2,length(subject),s)
%     [~,ch] = megPlotMap(sl_signal,climsSL,gcf,'bipolar');
%     hold on; contour(c2.XData,c2.YData, c2.ZData,3, 'k-'); colormap(bipolar);
%     set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);
    
%     subplot(2,length(subject),s+length(subject))
%     [~,ch] = megPlotMap(bb_signal,climsBB,gcf,'bipolar');
%     hold on; contour(c1.XData,c1.YData, c1.ZData,3, 'k-'); colormap(bipolar);
%     set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);
%     
    
    %% Calculate COD
    
%     absnorm = @(x, y, idx) abs(x) ./ norm([abs(x(idx)), abs(y(idx))]);
    normP = @(x, y, idx) (x) ./ norm([x(idx), y(idx)]);

    calcod = @(x, y, idx) 1 - sum((y(idx)-x(idx)).^2)./(sum(x(idx).^2));
    
    % Grab subject's data
    idx_d = isfinite(sl_signal);
    assert(isequal(idx_d, isfinite(bb_signal)));
    
    idx_p = isfinite(w.V1c(s,1:157));
    assert(isequal(idx_p, isfinite(w.V1i(s,1:157))));

    % Normalize data
    allDataSL_norm(s,:) = sl_signal./norm(sl_signal(idx_d));
    allDataBB_norm(s,:) = bb_signal./norm(bb_signal(idx_d));
    
    allPredictionCoherent_norm(s,:) = normP(w.V1c(s,1:157), w.V1i(s,1:157), idx_p);
    allPredictionIncoherent_norm(s,:) = normP(w.V1i(s,1:157), w.V1c(s,1:157), idx_p);

    
end


%% SAVING

if saveFigures % use different function to save figures, since figurewrite crashes with many subplots containing many data points
    
    set(0, 'currentfigure', fH1);
%     figurewrite(fullfile(figureDir,'Figure3A_predictionV123VsDataIndividuals_matched'),[],0,'.',1);
    hgexport(fH1, fullfile(figureDir, 'Figure3A_predictionV123VsDataIndividuals_matched.eps'))
    set(0, 'currentfigure', fH3);
%     figurewrite(fullfile(figureDir,'Figure3B_predictionV123VsDataIndividuals_matched'),[],0,'.',1);
    hgexport(fH3, fullfile(figureDir, 'Figure3B_predictionV123VsDataIndividuals_unmatched.eps'))
    
end

close all;

%% Calculate CoD
    
for d = 1:length(subject)
    
    thisDataSL = allDataSL_norm(d,:);
    thisDataBB = allDataBB_norm(d,:);
    
    idx = isfinite(thisDataSL);

    
    for p = 1:length(subject)
        
        thisPredictionCoherent = allPredictionCoherent_norm(p,:);
        thisPredictionIncoherent = allPredictionIncoherent_norm(p,:);
                        
        codSLCoherent(d, p, :) =  calcod(thisDataSL, thisPredictionCoherent, idx);
        codSLIncoherent(d, p, :) =  calcod(thisDataSL, thisPredictionIncoherent, idx);
        codBBCoherent(d, p, :) =  calcod(thisDataBB, thisPredictionCoherent, idx);
        codBBIncoherent(d, p, :) =  calcod(thisDataBB, thisPredictionIncoherent, idx);

        
    end
    
end


%%

slDataCoherentPred_subjectsMatched = diag(codSLCoherent);
slDataCoherentPred_subjectsNotMatched = mean(reshape(codSLCoherent(~eye(length(subject))),[length(subject)-1,length(subject)]))';

% bbDataIncoherentPred_subjectMatched = diag(codBBIncoherent);
% bbDataIncoherentPred_subjectNotMatched = mean(reshape(codBBIncoherent(~eye(6)),[5,6]))';

slDataIncoherentPred_subjectsMatched = diag(codSLIncoherent);
% bbDataCoherentPred_subjectsMatched = diag(codBBCoherent);


allData = cat(2, slDataCoherentPred_subjectsMatched, ...
                 slDataIncoherentPred_subjectsMatched, ...
                 slDataCoherentPred_subjectsNotMatched);                 
%                  bbDataIncoherentPred_subjectMatched, ...
%                  bbDataIncoherentPred_subjectNotMatched, ...                 
%                  bbDataCoherentPred_subjectsMatched);
                 
labels = {'Same subject - SL & Uniform', ...
          'Same subject - SL & Random', ...
          'Diff subject - SL & Uniform'};
%       , ... 
%           'Same subject - BB & Random', ...
%           'Diff subject - BB & Random', ...          
%            ...
%           'Same subject - BB & Uniform'};
          
colors = [0 0 0];

figure; set(gcf, 'Color', 'w', 'Position', [ 1000, 781, 335, 557])
boxplot(allData, 'Colors', colors, 'BoxStyle','outline', 'Widths',0.2, 'MedianStyle','line'); hold on
set(gca, 'YLim', [0 1]);
box off; set(gca, 'TickDir', 'out', 'TickLength',[0.015 0.015],'FontSize',20)
ylabel('Coefficient of Determination','FontSize',20); 
set(gca,'XTickLabel',labels); set(gca,'XTickLabelRotation',45);

hgexport(gcf, fullfile(figureDir, 'Figure3C_CoD_SL_UniformPrediction.eps'))








