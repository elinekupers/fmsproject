function makeFigure4


% % Add Fieldtrip and plotting function if not added yet.
% if ~exist('megPlotMap','file')
%     addpath(genpath('~/matlab/git/denoiseproject'))
% end
%
% if ~exist('ft_prepare_layout','file')
%     tbUse('Fieldtrip')
% end

% Set up paths
figureDir       = fullfile(fmsRootPath, 'figures'); % Where to save images?
saveFigures     = true;     % Save figures in the figure folder?

fH1 = figure(1); clf; set(fH1, 'Position', [1 1 1600 800], 'Name','Figure 4, Data against model predictions V1');

rg = [-1 1]*10E-5;

dataDir         = fullfile('~/matlab/git/denoiseproject/', 'analysis', 'data');    % Where to save data?


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
    
    d = dir(fullfile(bs_db, project_name, 'data', subject{s}));
    if strcmp(subject{s},'wl_subj002')
        data_dir = fullfile(bs_db, project_name, 'data', subject{s}, d(end-1).name);
    else
        data_dir = fullfile(bs_db, project_name, 'data', subject{s}, d(end).name);
    end
    
    anat_dir = fullfile(bs_db, project_name, 'anat', subject{s});
    
    %% 1. Load relevant matrices
    % Load Gain matrix created by brainstorm
    headmodel = load(fullfile(data_dir, 'headmodel_surf_os_meg.mat'));
    G = headmodel.Gain; % Gain matrix, [Nsensors x 3*Nvertices]
    G_constrained = bst_gain_orient(G, headmodel.GridOrient); % Contrained gain matrix [Nsensors x Nsources], equivalent to size pial cortex [1x15002]
    G_constrained = (G_constrained);
    
    % Load V1-3 template in downsampled brainstorm format (computed by interp_retinotopy.m)
    template = load(fullfile(anat_dir, 'areas_overlay.mat')); % [1xNsources] Every value between [-3,3] is inside V1-3, zeros refer to outside of visual cortex
    %     V123template = abs(template.sub_bs_areas)>0;
    V1template = abs(template.sub_bs_areas)==1;
    %     V2template = abs(template.sub_bs_areas)==2;
    %     V3template = abs(template.sub_bs_areas)==3;
    
    V123PhaseScrambledTemplate = load(fullfile(anat_dir, sprintf('areas_overlay_%s.mat',iterations))); % [1xNsources] Every value between [-3,3] is inside V1-3, zeros refer to outside of visual cortex
    
    % Load V1-3 eccentricity in downsampled brainstorm format
    eccen    = load(fullfile(anat_dir, 'eccen_overlay.mat')); % [1xNsources] Every value from [1 3] is inside V1-3, zeros refer to outside of visual cortex
    polarang = load(fullfile(anat_dir, 'angle_overlay.mat'));
    polarang.sub_bs_angle_rad = pi/180*(90-polarang.sub_bs_angle);
    
    % Restrict V123 to 11 degrees eccentricity (stimulus diameter = 22)
    %     V123templateStimEccen = V123template.*(eccen.sub_bs_eccen<=11);
    V1templateStimEccen = V1template.*(eccen.sub_bs_eccen<=11);
    %     V2templateStimEccen = V2template.*(eccen.sub_bs_eccen<=11);
    %     V3templateStimEccen = V3template.*(eccen.sub_bs_eccen<=11);
    
    %     V123PhaseScrambledTemplateStimEccen = V123PhaseScrambledTemplate.sub_bs_areas.*(eccen.sub_bs_eccen<=11);
    V1PhaseScrambledTemplateStimEccen   = V123PhaseScrambledTemplate.sub_bs_areas.*V1templateStimEccen;
    %     V2PhaseScrambledTemplateStimEccen   = V123PhaseScrambledTemplate.sub_bs_areas.*V2templateStimEccen;
    %     V3PhaseScrambledTemplateStimEccen   = V123PhaseScrambledTemplate.sub_bs_areas.*V3templateStimEccen;
    
    
    %% 2. Compute forward solution (Method 1)
    % Compute the sensor weights, w, from V1-V3 using the contrained gain field (forward model)
    %     w = G_constrained*V1template'; %  Nsensors x 1;
    
    phAmp2complex = @(r,th) r .* exp(1i*th);
    % Make a complex number with amplitudes and scrambled phases:
    V1complex = phAmp2complex(V1templateStimEccen,V1PhaseScrambledTemplateStimEccen);
    
    wComplexV1 = G_constrained*V1complex';
    
    % Compute sensor weights for restricted by stimulus eccentricity
    w1_stimsize = G_constrained*V1templateStimEccen'; %  Nsensors x 1;
    
    

    
    %% Data
    
    % Get stimulus locked and broadband response
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
    
    % Load denoised data of example subject
    [data] = prepareData(dataDir,whichSubject,5);
    bb = data{1};
    sl = data{2};
    
    % get stimulus-locked snr
    sl_signal = getsignalnoise(sl.results.origmodel(1),contrasts, 'SNR',sl.badChannels);
    % get broadband snr for before and after denoising
    bb_signal = getsignalnoise(bb.results.origmodel(1),  contrasts, 'SNR',bb.badChannels);
    
    
    sl_signal = to157chan(sl_signal,~sl.badChannels,'nans');
    bb_signal = to157chan(bb_signal,~bb.badChannels,'nans');
    
    
    %% Plotting
    figure(1); clf;
    ax1 =subplot(211);
    megPlotMap(abs(w1_stimsize(1:157)),rg,[],bipolar,[],[],[],'isolines', 3);
    c1 = findobj(ax1.Children,'Type','Contour');
    
    ax2 =subplot(212);
    megPlotMap(abs(wComplexV1(1:157)),0.5*rg,[],bipolar,[],[],[],'isolines', 3);
    c2 = findobj(ax2.Children,'Type','Contour');
    
    
    
    figure(2);
    subplot(2,length(subject),s)
    [~,ch] = megPlotMap(sl_signal,climsSL,gcf,'bipolar');
    hold on; contour(c1.XData,c1.YData, c1.ZData,3, 'k-'); colormap(bipolar); 
    set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);
    
    subplot(2,length(subject),s+length(subject))
    [~,ch] = megPlotMap(bb_signal,climsBB,gcf,'bipolar');
    hold on; contour(c2.XData,c2.YData, c2.ZData,3, 'k-'); colormap(bipolar); 
    set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);
    
end

%% SAVING

if saveFigures % use different function to save figures, since figurewrite crashes with many subplots containing many data points
%     figurewrite(fullfile(figureDir, ['predictionV1VsDataIndividuals']), [],0,'.',1);
end
