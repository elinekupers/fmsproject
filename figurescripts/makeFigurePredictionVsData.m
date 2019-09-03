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
dataDir         = fullfile(fmsRootPath, 'data');    % Where to get data?
projectName     = 'SSMEG';                          % Define Brainstorm project name, for subject and data/anatomy folders
saveFigures     = true;                             % Save figures in the figure folder?

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

% What type of data to use? 
dataType        = 'amplitudes'; % can be 'SNR' or 'amplitudes'
area            = 'V123'; % can be 'V123', 'V1', 'V2', 'V3' 
eccenLimitDeg   = [0.18 11]; % deg 

% What's the plotting range
climsSL         = [-40,40];
climsBB         = [-6,6];
climsContour    = [-1 1]*1E-4;

% Number of iterations for the random coherence prediction of the forward model
n        	= 10;        % number of timepoints (ms)
nrEpochs    = 1000;      % number of epochs
theta       = 0;         % von mises mean, equal for three distributions (syn, asyn and mix)
kappa.syn   = 100*pi;    % very narrow von Mises
kappa.asyn  = 0;         % very broad (uniform) von Mises
kappa.mix   = 0.27*pi;   % in-between width size von Mises (note: we are not plotting these values for this figure)

% Define vector that can truncate number of sensors 
keep_sensors    = logical([ones(157,1); zeros(192-157,1)]); % NB: Figure out a more generic way to define keep_sensors

% Predefine figures
fH1 = figure(1); clf; set(fH1, 'Position', [1 1 2550 1300], 'Name','Figure XX, Data against model predictions V1 - matched');
nrows = 4; 
ncols = 6;

% Loop over subjects to get predictions and data
for s = 1:length(subject)
    
    d = dir(fullfile(bsDB, projectName, 'data', subject{s}, 'R*'));
    BSdataDir = fullfile(d(1).folder, d(1).name);    
    BSanatDir = fullfile(bsDB, projectName, 'anat', subject{s});
    
    %% 1. Load relevant matrices
    % Get gain matrix
    G_constrained = getGainMatrix(BSdataDir, keep_sensors);

    % Get V1 template limited to 11 degrees eccentricity
    template = getTemplate(BSanatDir, area, eccenLimitDeg);

    % Simulate coherent, in between or mixture, adn incoherent source time 
    % series and compute predictions from forward model (w)
    tmp = getForwardModelPredictions(G_constrained, template.([area '_StimEccen']), [], n, nrEpochs, theta, kappa);
      
    % Compute amplitude across time
    amps.c = abs(fft(tmp.c,[],2));
    amps.i = abs(fft(tmp.i,[],2));
    amps.m = abs(fft(tmp.m,[],2));
    
    % Compute mean weights across epochs at input frequency
    w.V1c(s,:) = mean(amps.c(:,2,:),3);
    w.V1i(s,:) = mean(amps.i(:,2,:),3);
    w.V1m(s,:) = mean(amps.m(:,2,:),3);
    
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
        sl(s,:) = data{1};
        bb(s,:) = data{2};
    else
        if strcmp(subject{s},'wlsubj059')
            sl(s,:) = nanmean(data.sl.full_coherent,1) - nanmean(data.sl.blank_coherent,1);
        else
            sl(s,:) = nanmean(data.sl.full,1) - nanmean(data.sl.blank,1);
        end
        bb(s,:) = nanmean(data.bb.full,1) - nanmean(data.bb.blank,1);
        
        % Check for range (fT versus T)
        if max(sl(s,:)) < 1^-14
            sl(s,:) = sl(s,:) .* 10^15; % from T -> fT
            bb(s,:) = bb(s,:) .* 10^15 .* 10^15; % from T^2 -> fT.^2
        end
        
        climsSL = [-1 1]*prctile(sl(s,:), 97.5);
        climsBB = [-1 1]*prctile(bb(s,:), 97.5);
 
    end
       
    %% 4. Plotting to get contour lines
    figure(2); clf;
    ax1 = subplot(211);
    megPlotMap(abs(w.V1c(s,:)),climsContour,[],bipolar,[],[],[],'isolines', 3);
    c1 = findobj(ax1.Children,'Type','Contour');
    
    ax2 = subplot(212);
    megPlotMap(abs(w.V1i(s,:)),0.5*climsContour,[],bipolar,[],[],[],'isolines', 3);
    c2 = findobj(ax2.Children,'Type','Contour');
    
    figure(1);
    ax3 = subplot(nrows,ncols,s);
    [~,ch] = megPlotMap(sl(s,:),climsSL,gcf,'bipolar');
    hold on; contour(c1.XData,c1.YData, c1.ZData,3, 'k-');  drawnow; 
    c3 = findobj(ax3.Children,'Type','Contour'); c3.LineWidth = 3;
    colormap(bipolar); title(sprintf('SL: S%d',s));
    set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);
    
    ax4 = subplot(nrows,ncols,s+length(subject));
    [~,ch] = megPlotMap(bb(s,:),climsBB,gcf,'bipolar');
    hold on; contour(c2.XData,c2.YData, c2.ZData,3, 'k-'); drawnow;
    c4 = findobj(ax4.Children,'Type','Contour'); c4.LineWidth = 3;
    colormap(bipolar); title(sprintf('BB: S%d',s));
    set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);   
    
    %% Calculate COD
    
%     absnorm = @(x, y, idx) abs(x) ./ norm([abs(x(idx)), abs(y(idx))]);
    normP = @(x, y, idx) (x) ./ norm([x(idx), y(idx)]);

    calcod = @(x, y, idx) 1 - sum((y(idx)-x(idx)).^2)./(sum(x(idx).^2));
    
    % Grab subject's data
    idx_d = isfinite(sl(s,:));
    assert(isequal(idx_d, isfinite(bb(s,:))));
    
    idx_p = isfinite(w.V1c(s,1:157));
    assert(isequal(idx_p, isfinite(w.V1i(s,1:157))));

    % Normalize data
    allDataSL_norm(s,:) = sl(s,:)./norm(sl(s,idx_d));
    allDataBB_norm(s,:) = bb(s,:)./norm(bb(s,idx_d));
    
    allPredictionCoherent_norm(s,:) = normP(w.V1c(s,1:157), w.V1i(s,1:157), idx_p);
    allPredictionIncoherent_norm(s,:) = normP(w.V1i(s,1:157), w.V1c(s,1:157), idx_p);

    
end

%% SAVE Figure

if saveFigures % use different function to save figures, since figurewrite crashes with many subplots containing many data points
    
    set(0, 'currentfigure', fH1);
%     figurewrite(fullfile(figureDir,'Figure3A_predictionV123VsDataIndividuals_matched'),[],0,'.',1);
    hgexport(fH1, fullfile(figureDir, sprintf('FigureXX_prediction_%s_%1.2f-%d_VsDataIndividuals_matched_%s.eps', area, eccenLimitDeg(1), eccenLimitDeg(2), dataType)))
    
end

close all;

%% Plot average across subjects
figure(4); clf;
ax1 = subplot(211);
megPlotMap(nanmean(abs(w.V1c),1),climsContour,[],bipolar,[],[],[],'isolines', 3);
c1 = findobj(ax1.Children,'Type','Contour');

ax2 =subplot(212);
megPlotMap(nanmean(abs(w.V1i),1),0.5*climsContour,[],bipolar,[],[],[],'isolines', 3);
c2 = findobj(ax2.Children,'Type','Contour');

figure(5);
subplot(211)
[~,ch] = megPlotMap(nanmean(sl,1),climsSL,gcf,'bipolar');
hold on; contour(c1.XData,c1.YData, c1.ZData,3, 'k-');  drawnow;
colormap(bipolar); title(sprintf('SL: N = %d',length(subject)));
set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);

subplot(212)
[~,ch] = megPlotMap(nanmean(bb,1),climsBB,gcf,'bipolar');
hold on; contour(c2.XData,c2.YData, c2.ZData,3, 'k-'); drawnow;
colormap(bipolar); title(sprintf('BB N = %d',length(subject)));
set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);  


mnDataSL_norm = nanmean(allDataSL_norm,1);
mnDataBB_norm = nanmean(allDataBB_norm,1);
mnPredictionCoherent = nanmean(allPredictionCoherent_norm,1);
mnPredictionIncoherent = nanmean(allPredictionIncoherent_norm,1);

if saveFigures % use different function to save figures, since figurewrite crashes with many subplots containing many data points    
    set(0, 'currentfigure', 5);
%     figurewrite(fullfile(figureDir,'Figure3A_predictionV123VsDataIndividuals_matched'),[],0,'.',1);
    hgexport(5, fullfile(figureDir, sprintf('Figure3A_prediction_%s_%1.2f-%d_VsDataAVERAGE_matched_%s.eps', area, eccenLimitDeg(1), eccenLimitDeg(2), dataType)))
    
end

close all;

%% Calculate CoD for individual subjects
    
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

%% for mean across subjects:
idx = isfinite(mnDataSL_norm);
mn_codSLCoherent   =  calcod(mnDataSL_norm, mnPredictionCoherent, idx);
mn_codSLIncoherent =  calcod(mnDataSL_norm, mnPredictionIncoherent, idx);
mn_codBBCoherent   =  calcod(mnDataBB_norm, mnPredictionCoherent, idx);
mn_codBBIncoherent =  calcod(mnDataBB_norm, mnPredictionIncoherent, idx);


%% Plot bar graph
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


fprintf('\nAVERAGE ACROSS SUBJECTS:\n')
fprintf('\nEFFECT OF SYNCHRONY when MATCHED:\n')
fprintf('CoD SL w/ coherent phase: %1.3f\n', mn_codSLCoherent)
fprintf('CoD BB w/ incoherent phase: %1.3f\n', mn_codBBIncoherent)

fprintf('\nAVERAGE ACROSS SUBJECTS:\n')
fprintf('\nEFFECT OF SYNCHRONY when NOT MATCHED:\n')
fprintf('CoD SL w/ incoherent phase: %1.3f\n', mn_codSLIncoherent)
fprintf('CoD BB w/ coherent phase: %1.3f\n', mn_codBBCoherent)



