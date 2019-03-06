function makeFigure_dataVsPrediction()

% This is a function to make a figure from the manuscript about forward
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

%% 0. Set up paths and define parameters
figureDir       = fullfile(fmsRootPath, 'figures'); % Where to save images?
saveFigures     = true;     % Save figures in the figure folder?

fH1 = figure(1); clf; set(fH1, 'Position', [1 1 1600 800], 'Name','Figure 4A, SL Data against coh model predictions V123');
fH3 = figure(3); clf; set(fH3, 'Position', [1 1 1600 800], 'Name','Figure 4B, BB Data against incoh model predictions V123');
fH4 = figure(4); clf; set(fH4, 'Position', [1 1 1600 800], 'Name','Figure 4c, SL Data against mix model predictions V123');

rg = [-1 1]*10E-4;

dataDir         = fullfile(fmsRootPath, 'data');    % Where to get data?

climsSL = [-20,20];
climsBB = [-6,6];

% Use coherent or incoherent spectrum to compute the SL signal
useSLIncohSpectrum = false;

% Define vector that can truncate number of sensors
keep_sensors = logical([ones(157,1); zeros(192-157,1)]); % Note: Figure out a more generic way to define keep_sensors

% What visual area to use?
area            = 'all'; % Choose between 'V1' or 'all' (=V1-V3);

% Number of iterations for the random coherence prediction of the forward
% model
n           = 10;         % number of timepoints (ms)
nrEpochs    = 1000;       % number of epochs
theta       = 0;          % von mises mean of three distributions
kappa.coh   = 10*pi;
kappa.incoh = 0;
kappa.mix   = 0.5*pi;

%% 0. Define paths
% Path to brainstorm database
bsDB = '/Volumes/server/Projects/MEG/brainstorm_db/';

% Define project name, subject and data/anatomy folders
projectName = 'SSMEG';

% Which subjects to average?
subject   = {'wlsubj002','wlsubj004','wlsubj005','wlsubj006','wlsubj010','wlsubj011','wlsubj048', 'wlsubj046','wlsubj039','wlsubj059', 'wlsubj067', 'wlsubj070'};

for s = 1:length(subject)
    
    d = dir(fullfile(bsDB, projectName, 'data', subject{s}, 'R*'));
    bsData = fullfile(d(1).folder, d(1).name);
    bsAnat = fullfile(bsDB, projectName, 'anat', subject{s});
    
    % Get constrained Gain matrix
    G_constrained = getGainMatrix(bsData, keep_sensors);
        
    % Get V1 template limited to 11 degrees eccentricity
    template = getTemplate(bsAnat, area, 11);

    % Simulate coherent, in between or mixture, adn incoherent source time
    % series and compute predictions from forward model (w)
    if strcmp(area, 'all')
        tmp = getForwardModelPredictions(G_constrained, template.V123StimEccen, [], n, nrEpochs, theta, kappa);
    else
        tmp = getForwardModelPredictions(G_constrained, template.V1StimEccen, [], n, nrEpochs, theta, kappa);
    end

    % Compute amplitude across time
    amps.c = abs(fft(tmp.c,[],2));
    amps.i = abs(fft(tmp.i,[],2));
    amps.m = abs(fft(tmp.m,[],2));

    % Compute mean weights across epochs at input frequency
    w.V1c(s,:) = mean(amps.c(:,2,:),3);
    w.V1i(s,:) = mean(amps.i(:,2,:),3);
    w.V1m(s,:) = mean(amps.m(:,2,:),3);
    
    
    %% Data
    
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
    

    %% Plotting
    figure(2); clf;
    ax1 =subplot(311);
    megPlotMap(w.V1c(s,:),rg,[],bipolar,[],[],[],'isolines', 3);
    c1 = findobj(ax1.Children,'Type','Contour');
    title('Kappa = 10*pi')
    
    ax2 =subplot(312);
    megPlotMap(w.V1i(s,:),rg,[],bipolar,[],[],[],'isolines', 3);
    c2 = findobj(ax2.Children,'Type','Contour');
    title('Kappa = 0')
    
    ax3 =subplot(313);
    megPlotMap(w.V1m(s,:),rg,[],bipolar,[],[],[],'isolines', 3);
    c3 = findobj(ax3.Children,'Type','Contour');  
    title('Kappa = 0.5*pi')
    
    hgexport(gcf, fullfile(figureDir, subject{s}, 'FigureXX_predictionV123_contours.eps'))

    
    figure(1);
    subplot(2,ceil(length(subject)/2),s)
    [~,ch] = megPlotMap(diffFullBlankSL(s,:),climsSL,gcf,'bipolar');
    hold on; contour(c1.XData,c1.YData, c1.ZData,3, 'k-'); colormap(bipolar);
    set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);
    title(sprintf('S%d',s))

    figure(3);
    subplot(2,ceil(length(subject)/2),s)
    [~,ch] = megPlotMap(diffFullBlankBB(s,:),climsBB,gcf,'bipolar');
    hold on; contour(c2.XData,c2.YData, c2.ZData,3, 'k-'); colormap(bipolar);
    set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);
    title(sprintf('S%d',s))
    
    figure(4);
    subplot(2,ceil(length(subject)/2),s)
    [~,ch] = megPlotMap(diffFullBlankSL(s,:),climsSL,gcf,'bipolar');
    hold on; contour(c3.XData,c3.YData, c3.ZData,3, 'k-'); colormap(bipolar);
    set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);
    title(sprintf('S%d',s))
    
%     figure(3);
%     subplot(2,length(subject),s)
%     [~,ch] = megPlotMap(diffFullBlankSL(s,:),climsSL,gcf,'bipolar');
%     hold on; contour(c2.XData,c2.YData, c2.ZData,3, 'k-'); colormap(bipolar);
%     set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);
%     
%     subplot(2,length(subject),s+length(subject))
%     [~,ch] = megPlotMap(diffFullBlankBB(s,:),climsBB,gcf,'bipolar');
%     hold on; contour(c1.XData,c1.YData, c1.ZData,3, 'k-'); colormap(bipolar);
%     set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);
%     
    
    %% Calculate COD
%     
%     absnorm = @(x, idx) abs(x) ./ norm(abs(x(idx)));
%     calcod = @(x, y, idx) 1 - sum((y(idx)-x(idx)).^2)./(sum(x(idx).^2));
%     
%     % Grab subject's data
%     idx_d = isfinite(diffFullBlankBB(s,:));
%     assert(isequal(idx_d, isfinite(diffFullBlankBB(s,:))));
%     
%     idx_p = isfinite(w.V1c(s,:));
%     assert(isequal(idx_p, isfinite(w.V1i(s,:))));
% 
%     
%     % Normalize data
%     allDataSL_norm(s,:) = diffFullBlankSL(s,:)./norm(diffFullBlankSL(s,idx_d));
%     allDataBB_norm(s,:) = diffFullBlankBB(s,:)./norm(diffFullBlankBB(s,idx_d));
%     
%     allPredictionUniform_norm(s,:) = absnorm(w.V1c(s,:), idx_p);
%     allPredictionRandom_norm(s,:) = absnorm(w.V1i(s,:), idx_p);

    
end


%% SAVING

if saveFigures % use different function to save figures, since figurewrite crashes with many subplots containing many data points
    hgexport(fH1, fullfile(figureDir, 'FigureXXA_predictionV123VsDataIndividuals_Coh_spectrumCoh.eps'))
    hgexport(fH3, fullfile(figureDir, 'FigureXXB_predictionV123VsDataIndividuals_Incoh.eps'))
    hgexport(fH4, fullfile(figureDir, 'FigureXXC_predictionV123VsDataIndividuals_Mix_spectrumCoh.eps'))

    %     figurewrite(fullfile(figureDir, ['predictionV1VsDataIndividuals']), [],0,'.',1);
end

close all;

%% Calculate CoD
    
for d = 1:6
    
    thisDataSL = allDataSL_norm(d,:);
    thisDataBB = allDataBB_norm(d,:);
    
    idx = isfinite(thisDataSL);

    
    for p = 1:6
        
        thisPredictionUniform = allPredictionUniform_norm(p,:);
        thisPredictionRandom = allPredictionRandom_norm(p,:);
                        
        codSLUniform(d, p, :) =  calcod(thisDataSL, thisPredictionUniform, idx);
        codSLRandom(d, p, :) =  calcod(thisDataSL, thisPredictionRandom, idx);
        codBBUniform(d, p, :) =  calcod(thisDataBB, thisPredictionUniform, idx);
        codBBRandom(d, p, :) =  calcod(thisDataBB, thisPredictionRandom, idx);

        
    end
    
end


% Plot CoD 
figure; imagesc(codSLUniform); colormap gray; axis square; xlabel('Data SL'); ylabel('Forward model V1 Uniform phase Prediction'); colorbar; makeprettyaxes(gca,9,9)
title(sprintf('R squared - mean diag: %1.2f, mean off diag: %1.2f',mean(diag(codSLUniform)), mean(codSLUniform(~eye(6)))));
range = get(gca,'CLim');

printnice(gcf, [1 300], fullfile(figureDir),'Figure4C_CoD_SL_UniformPrediction');
figurewrite(fullfile(figureDir, 'Figure4C_CoD_SL_UniformPrediction.eps'))

figure; imagesc(codSLRandom); colormap gray; axis square; xlabel('Data SL'); ylabel('Forward model V1 Randorm phase Prediction'); colorbar; makeprettyaxes(gca,9,9)
title(sprintf('R squared - mean diag: %1.2f, mean off diag: %1.2f',mean(diag(codSLRandom)), mean(codSLRandom(~eye(6)))))
set(gca, 'CLim', range)

printnice(gcf, [1 300], fullfile(figureDir),'Figure4D_CoD_SL_ScrambledPrediction');
hgexport(gcf, fullfile(figureDir, 'Figure4D_CoD_SL_ScrambledPrediction.eps'))

%%
figure;
plot(0.5*ones(length(diag(codSLUniform)),1),diag(codSLUniform), 'ko','MarkerSize',15); hold on
plot([0.35 0.65],[mean(diag(codSLUniform)), mean(diag(codSLUniform))],'r', 'LineWidth',4);

plot(ones(length(codSLUniform(~eye(6))),1),codSLUniform(~eye(6)), 'ko','MarkerSize',15); hold on
plot([.85 1.15],[mean(codSLUniform(~eye(6))), mean(codSLUniform(~eye(6)))],'r', 'LineWidth',4);

xlim([0 1.5]); ylim([0 1]); box off; set(gca, 'TickDir', 'out', 'XTick',[0.5 1], 'XTickLabel', {'on','off'}, 'TickLength',[0.015 0.015],'FontSize',20)
ylabel('Coefficient of Determination','FontSize',20); title('SL data and Uniform prediction')

hgexport(gcf, fullfile(figureDir, 'Figure4E_CoD_SL_UniformPrediction.eps'))


figure;
plot(0.5*ones(length(diag(codSLRandom)),1),diag(codSLRandom), 'ko','MarkerSize',15); hold on
plot([0.35 0.65],[mean(diag(codSLRandom)), mean(diag(codSLRandom))],'r', 'LineWidth',4);

plot(ones(length(codSLRandom(~eye(6))),1),codSLRandom(~eye(6)), 'ko','MarkerSize',15); hold on
plot([.85 1.15],[mean(codSLRandom(~eye(6))), mean(codSLRandom(~eye(6)))],'r', 'LineWidth',4);

xlim([0 1.5]); ylim([0 1]); box off; set(gca, 'TickDir', 'out', 'XTick',[0.5 1], 'XTickLabel', {'on','off'}, 'TickLength',[0.015 0.015],'FontSize',20)
ylabel('Coefficient of Determination','FontSize',20); title('SL data and Random prediction')

hgexport(gcf, fullfile(figureDir, 'Figure4F_CoD_SL_RandomPrediction.eps'))












