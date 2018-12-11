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
projectName     = 'SSMEG';

% Which subjects to average?
%   Full  only: 'wlsubj048', 'wlsubj046','wlsubj039','wlsubj059', 'wlsubj067'
%   Full, Left, Right: 'wlsubj002','wlsubj004','wlsubj005','wlsubj006','wlsubj010','wlsubj011'
subject         = {'wlsubj002','wlsubj004','wlsubj005','wlsubj006','wlsubj010','wlsubj011','wlsubj048', 'wlsubj046','wlsubj039','wlsubj059', 'wlsubj067', 'wlsubj070'};

% What type of data to use? 
dataType        = 'amplitudes'; % can be 'SNR' or 'amplitudes'
area            = 'all'; % can be 'all' (= V1-V3), or 'V1' 

% What's the plotting range
climsSL         = [-20,20];
climsBB         = [-6,6];
climsPred       = [-1 1]*1E-4;

% Number of iterations for the random coherence prediction of the forward model
n               = 10;     % number of timepoints (ms)
nrEpochs        = 1;      % number of epochs

% Define vector that can truncate number of sensors 
keep_sensors    = logical([ones(157,1); zeros(192-157,1)]); % Note: Figure out a more generic way to define keep_sensors

% Predefine figures
fH1 = figure(1); clf; set(fH1, 'Position', [1 1 2550 1300], 'Name','Figure 3A, Data against model predictions V1 - matched');
if length(subject) <= 6; nrows = 2; ncols = 6; else; nrows = 4; ncols = 6; end

% Loop over subjects to get predictions and data
for s = 1:length(subject)
    
    d = dir(fullfile(bsDB, projectName, 'data', subject{s}, 'R*'));
    BSdataDir = fullfile(d(1).folder, d(1).name);    
    BSanatDir = fullfile(bsDB, projectName, 'anat', subject{s});
    
    %% 1. Load relevant matrices
    % Get gain matrix
    G_constrained = getGainMatrix(BSdataDir, keep_sensors);

    % Get V1 template limited to 11 degrees eccentricity
    template = getTemplate(BSanatDir, area, 11);

    % Simulate coherent and incoherent source time series and compute
    % predictions from forward model (w)
    if strcmp(area, 'all')
        tmp = getForwardModelPredictions(G_constrained, template.V123StimEccen, [], n, nrEpochs);
    else
       tmp = getForwardModelPredictions(G_constrained, template.V1StimEccen, [], n, nrEpochs);
    end
    
    % Take mean across epochs
    w.V1c(s,:) = mean(abs(tmp.c),2);
    w.V1i(s,:) = mean(abs(tmp.i),2);
    
    %% 3. Data
    
    % Get stimulus locked and broadband response
    switch subject{s}
        % Go from subject to session nr
        case 'wlsubj002'
            whichSession = 2;
        case 'wlsubj004'
            whichSession = 7;
        case 'wlsubj005'
            whichSession = 8;
        case 'wlsubj006'
            whichSession = 1;
        case 'wlsubj010'
            whichSession = 6;
        case 'wlsubj011'
            whichSession = 5;
        case 'wlsubj048'
            whichSession = 9; % Full field Only
        case 'wlsubj046'
            whichSession = 10; % Full field Only
        case 'wlsubj039'
            whichSession = 11; % Full field Only
        case 'wlsubj059'
            whichSession = 12; % Full field Only
        case 'wlsubj067'
            whichSession = 13; % Full field Only
        case 'wlsubj070'
            whichSession = 14; % Full field Only
    end
    
    % Load denoised data of example subject
    data = loadData(fullfile(dataDir, subject{s}),whichSession, dataType);
    if strcmp(dataType, 'SNR')
        sl = data{1};
        bb = data{2};
    else
        sl = nanmean(data.sl.full_coherent,1) - nanmean(data.sl.blank_coherent,1);
        bb = nanmean(data.bb.full,1) - nanmean(data.bb.blank,1);
        climsSL = [-1 1]*prctile(sl, 97.5);
        climsBB = [-1 1]*prctile(bb, 97.5);
    end
       
    %% 4. Plotting to get contour lines
    figure(2); clf;
    ax1 = subplot(211);
    megPlotMap(abs(w.V1c(s,:)),climsPred,[],bipolar,[],[],[],'isolines', 3);
    c1 = findobj(ax1.Children,'Type','Contour');
    
    ax2 =subplot(212);
    megPlotMap(abs(w.V1i(s,:)),0.5*climsPred,[],bipolar,[],[],[],'isolines', 3);
    c2 = findobj(ax2.Children,'Type','Contour');
    
    figure(1);
    subplot(nrows,ncols,s)
    [~,ch] = megPlotMap(sl,climsSL,gcf,'bipolar');
    hold on; contour(c1.XData,c1.YData, c1.ZData,3, 'k-'); 
    colormap(bipolar); title(sprintf('SL: S%d',s));
    set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);
    
    subplot(nrows,ncols,s+length(subject))
    [~,ch] = megPlotMap(bb,climsBB,gcf,'bipolar');
    hold on; contour(c2.XData,c2.YData, c2.ZData,3, 'k-');
    colormap(bipolar); title(sprintf('BB: S%d',s));
    set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);   
    
    %% Calculate COD
    
%     absnorm = @(x, y, idx) abs(x) ./ norm([abs(x(idx)), abs(y(idx))]);
    normP = @(x, y, idx) (x) ./ norm([x(idx), y(idx)]);

    calcod = @(x, y, idx) 1 - sum((y(idx)-x(idx)).^2)./(sum(x(idx).^2));
    
    % Grab subject's data
    idx_d = isfinite(sl);
    assert(isequal(idx_d, isfinite(bb)));
    
    idx_p = isfinite(w.V1c(s,1:157));
    assert(isequal(idx_p, isfinite(w.V1i(s,1:157))));

    % Normalize data
    allDataSL_norm(s,:) = sl./norm(sl(idx_d));
    allDataBB_norm(s,:) = bb./norm(bb(idx_d));
    
    allPredictionCoherent_norm(s,:) = normP(w.V1c(s,1:157), w.V1i(s,1:157), idx_p);
    allPredictionIncoherent_norm(s,:) = normP(w.V1i(s,1:157), w.V1c(s,1:157), idx_p);

    
end

%% SAVING

if saveFigures % use different function to save figures, since figurewrite crashes with many subplots containing many data points
    
    set(0, 'currentfigure', fH1);
%     figurewrite(fullfile(figureDir,'Figure3A_predictionV123VsDataIndividuals_matched'),[],0,'.',1);
    hgexport(fH1, fullfile(figureDir, sprintf('Figure3A_prediction_%s_VsDataIndividuals_matched_%s.eps', area, dataType)))
    
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

% Plot bar graph
slDataCoherentPred_subjectsMatched = diag(codSLCoherent);
bbDataInCoherentPred_subjectsMatched = diag(codBBIncoherent);

slDataCoherentPred_subjectsNotMatched = mean(reshape(codSLCoherent(~eye(length(subject))),[length(subject)-1,length(subject)]))';
slDataIncoherentPred_subjectsMatched = diag(codSLIncoherent);

allData = cat(2, slDataCoherentPred_subjectsMatched, ...
                 bbDataInCoherentPred_subjectsMatched, ...
                 slDataIncoherentPred_subjectsMatched, ...
                 slDataCoherentPred_subjectsNotMatched);                 
                 
labels = {'Same subject - SL & Uniform', ...
          'Same subject - BB & random', ...
          'Same subject - SL & Random', ...
          'Diff subject - SL & Uniform'};
colors = [0 0 0];

figure; set(gcf, 'Color', 'w', 'Position', [ 1000, 781, 335, 557])
boxplot(allData, 'Colors', colors, 'BoxStyle','outline', 'Widths',0.2, 'MedianStyle','line'); hold on
set(gca, 'YLim', [-1.5 1]);
box off; set(gca, 'TickDir', 'out', 'TickLength',[0.015 0.015],'FontSize',20)
ylabel('Coefficient of Determination','FontSize',20); 
set(gca,'XTickLabel',labels); set(gca,'XTickLabelRotation',45);

hgexport(gcf, fullfile(figureDir, sprintf('Figure3C_CoD_SL_UniformPrediction_%s_%s.eps', area, dataType)))

fprintf('\nEFFECT OF SYNCHRONY when MATCHED:\n')
fprintf('CoD SL w/ coherent phase: Median: %1.3f, Mean (+/- se) = %1.3f (+/- %1.3f)\n', median(slDataCoherentPred_subjectsMatched), mean(slDataCoherentPred_subjectsMatched), (std(slDataCoherentPred_subjectsMatched))./sqrt(length(slDataCoherentPred_subjectsMatched)))
fprintf('CoD BB w/ incoherent phase: Median: %1.3f, Mean (+/- se) = %1.3f (+/- %1.3f)\n', median(bbDataInCoherentPred_subjectsMatched), mean(bbDataInCoherentPred_subjectsMatched), (std(bbDataInCoherentPred_subjectsMatched))./sqrt(length(bbDataInCoherentPred_subjectsMatched)))

fprintf('\nEFFECT OF SYNCHRONY when NOT MATCHED:\n')
fprintf('CoD SL w/ incoherent phase: Median: %1.3f, Mean (+/- se) = %1.3f (+/- %1.3f)\n', median(slDataIncoherentPred_subjectsMatched), mean(slDataIncoherentPred_subjectsMatched), (std(slDataIncoherentPred_subjectsMatched))./sqrt(length(slDataIncoherentPred_subjectsMatched)))


fprintf('\nEFFECT OF ANATOMY when NOT MATCHED:\n')
fprintf('CoD SL w/ coherent phase, but wrong anatomy: Median: %1.3f, Mean (+/- se) = %1.3f (+/- %1.3f)\n', median(slDataCoherentPred_subjectsNotMatched), mean(slDataCoherentPred_subjectsNotMatched), (std(slDataCoherentPred_subjectsNotMatched))./sqrt(length(slDataCoherentPred_subjectsNotMatched)))





