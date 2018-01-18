function makeSupplementaryFigure1()

% This is a function to make Supplementary Figure 1 from the manuscript about forward
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

% Set up paths
figureDir       = fullfile(fmsRootPath, 'figures'); % Where to save images?
saveFigures     = true;     % Save figures in the figure folder?
fH1 = figure(1); clf; set(fH1, 'Position', [1 1 1600 800], 'Name','Supl Figure 1A, Data against model predictions V1 - matched');
fH3 = figure(3); clf; set(fH3, 'Position', [1 1 1600 800], 'Name','Supl Figure 1B, Data against model predictions V1 - unmatched');

rg = [-1 1]*10E-5;

meg_data_dir         = fullfile(fmsRootPath, 'data');    % Where to get data?

contrasts = [1 0 0];
contrasts = bsxfun(@rdivide, contrasts, sqrt(sum(contrasts.^2,2)));
climsSL = [-20,20];
climsBB = [-6,6];

%% 0. Define paths
% Path to brainstorm database
bs_db = '/Volumes/server/Projects/MEG/brainstorm_db/';

% Define project name, subject and data/anatomy folders
project_name = 'SSMEG';

% How many iterations of smoothing??
iterations = 'phase_0'; % pick between 'phase_1', 'phase_5', 'phase_10', 'phase_20', 'phase_nosmoothing', 'phase_1_normalized' where the number stands for the number of smoothing iterations

% Which subjects to average?
subject = {'wl_subj002','wl_subj004','wl_subj005','wl_subj006','wl_subj010','wl_subj011'};

for s = 1:length(subject)
    
    % find data directories
    d = dir(fullfile(bs_db, project_name, 'data', subject, 'R*'));

    bs_data_dir = fullfile(d(1).folder, d(1).name);
    anat_dir = fullfile(bs_db, project_name, 'anat', subject{s});
    
    %% 1. Load relevant matrices
    % Load Gain matrix created by brainstorm
    headmodel = load(fullfile(bs_data_dir, 'headmodel_surf_os_meg.mat'));
    G = headmodel.Gain; % Gain matrix, [Nsensors x 3*Nvertices]
    G_constrained = bst_gain_orient(G, headmodel.GridOrient); % Contrained gain matrix [Nsensors x Nsources], equivalent to size pial cortex [1x15002]
    G_constrained = (G_constrained);
    
    % Load V1-3 template in downsampled brainstorm format (computed by interp_retinotopy.m)
    template = load(fullfile(anat_dir, 'areas_overlay.mat')); % [1xNsources] Every value between [-3,3] is inside V1-3, zeros refer to outside of visual cortex
    V123template = abs(template.sub_bs_areas)>0;

    V123PhaseScrambledTemplate = load(fullfile(anat_dir, sprintf('areas_overlay_%s.mat',iterations))); % [1xNsources] Every value between [-3,3] is inside V1-3, zeros refer to outside of visual cortex
    
    % Load V1-3 eccentricity in downsampled brainstorm format
    eccen    = load(fullfile(anat_dir, 'eccen_overlay.mat')); % [1xNsources] Every value from [1 3] is inside V1-3, zeros refer to outside of visual cortex
    polarang = load(fullfile(anat_dir, 'angle_overlay.mat'));
    polarang.sub_bs_angle_rad = pi/180*(90-polarang.sub_bs_angle);
    
    % Restrict V123 to 11 degrees eccentricity (stimulus diameter = 22)
    V123templateStimEccen = V123template.*(eccen.sub_bs_eccen<=11);
    
    V123PhaseScrambledTemplateStimEccen = V123PhaseScrambledTemplate.sub_bs_areas.*(eccen.sub_bs_eccen<=11);
    
    
    %% 2. Compute forward solution (Method 1)
    % Compute the sensor weights, w, from V1-V3 using the contrained gain field (forward model)
    %     w = G_constrained*V1template'; %  Nsensors x 1;
    
    phAmp2complex = @(r,th) r .* exp(1i*th);
    % Make a complex number with amplitudes and scrambled phases:
    V123complex = phAmp2complex(V123templateStimEccen,V123PhaseScrambledTemplateStimEccen);
    
    wComplexV123 = G_constrained*V123complex';
    
    % Compute sensor weights for restricted by stimulus eccentricity
    w123_stimsize = G_constrained*V123templateStimEccen'; %  Nsensors x 1;
    

    %% Data
    
    % Get stimulus locked and broadband response
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
    
    % Load denoised data of example subject
    data = loadData(meg_data_dir,whichSubject);
    bb = data{1};
    sl = data{2};
    
    % get stimulus-locked snr
    sl_signal = getsignalnoise(sl.results.origmodel(1),contrasts, 'SNR',sl.badChannels);
    % get broadband snr for before and after denoising
    bb_signal = getsignalnoise(bb.results.origmodel(1),  contrasts, 'SNR',bb.badChannels);
    
    
    sl_signal = to157chan(sl_signal,~sl.badChannels,'nans');
    bb_signal = to157chan(bb_signal,~bb.badChannels,'nans');
    
    
    %% Plotting
    figure(2); clf;
    ax1 =subplot(211);
    megPlotMap(abs(w123_stimsize(1:157)),rg,[],bipolar,[],[],[],'isolines', 3);
    c1 = findobj(ax1.Children,'Type','Contour');
    
    ax2 =subplot(212);
    megPlotMap(abs(wComplexV123(1:157)),0.5*rg,[],bipolar,[],[],[],'isolines', 3);
    c2 = findobj(ax2.Children,'Type','Contour');
    
    
    
    figure(1);
    subplot(2,length(subject),s)
    [~,ch] = megPlotMap(sl_signal,climsSL,gcf,'bipolar');
    hold on; contour(c1.XData,c1.YData, c1.ZData,3, 'k-'); colormap(bipolar);
    set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);
    
    subplot(2,length(subject),s+length(subject))
    [~,ch] = megPlotMap(bb_signal,climsBB,gcf,'bipolar');
    hold on; contour(c2.XData,c2.YData, c2.ZData,3, 'k-'); colormap(bipolar);
    set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);
    
    figure(3);
    subplot(2,length(subject),s)
    [~,ch] = megPlotMap(sl_signal,climsSL,gcf,'bipolar');
    hold on; contour(c2.XData,c2.YData, c2.ZData,3, 'k-'); colormap(bipolar);
    set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);
    
    subplot(2,length(subject),s+length(subject))
    [~,ch] = megPlotMap(bb_signal,climsBB,gcf,'bipolar');
    hold on; contour(c1.XData,c1.YData, c1.ZData,3, 'k-'); colormap(bipolar);
    set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);
    
    
    %% Calculate COD
    
    absnorm = @(x, idx) abs(x) ./ norm(abs(x(idx)));
    calcod = @(x, y, idx) 1 - sum((y(idx)-x(idx)).^2)./(sum(x(idx).^2));
    
    % Grab subject's data
    idx_d = isfinite(sl_signal);
    assert(isequal(idx_d, isfinite(bb_signal)));
    
    idx_p = isfinite(w123_stimsize(1:157));
    assert(isequal(idx_p, isfinite(wComplexV123(1:157))));

    
    % Normalize data
    allDataSL_norm(s,:) = sl_signal./norm(sl_signal(idx_d));
    allDataBB_norm(s,:) = bb_signal./norm(bb_signal(idx_d));
    
    allPredictionUniform_norm(s,:) = absnorm(w123_stimsize(1:157), idx_p);
    allPredictionRandom_norm(s,:) = absnorm(wComplexV123(1:157), idx_p);

    
end


%% SAVING

if saveFigures % use different function to save figures, since figurewrite crashes with many subplots containing many data points
    hgexport(fH1, fullfile(figureDir, 'SupFigure1A_predictionV123VsDataIndividuals_matched.eps'))
    hgexport(fH3, fullfile(figureDir, 'SupFigure1B_predictionV123VsDataIndividuals_unmatched.eps'))
    
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


%%

slDataPredMatched = diag(codSLUniform);
slDataPredNotMatched = mean(reshape(codSLUniform(~eye(6)),[5,6]))';
bbDataPredMatched = diag(codSLRandom);

allData = cat(2, slDataPredMatched, bbDataPredMatched, slDataPredNotMatched);
labels = {'Matched Coherent Prediction', ...
          'Matched Incoherent Prediction', ...
          'Unmatched Coherent Prediction'};
          
colors = [0 0 0];

figure;
boxplot(allData, 'Colors', colors, 'BoxStyle','outline', 'Widths',0.2); hold on

ylim([0 1]); box off; set(gca, 'TickDir', 'out', 'XTick',[1:3], 'XTickLabel', {'on','off'}, 'TickLength',[0.015 0.015],'FontSize',20)
ylabel('Coefficient of Determination','FontSize',20); 
set(gca,'XTickLabel',labels); set(gca,'XTickLabelRotation',45);

hgexport(gcf, fullfile(figureDir, 'SupFigure1C_CoD_SL_UniformPrediction.eps'))








