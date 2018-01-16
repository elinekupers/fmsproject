function makeFigure5()

% This is a function to make Figure 5 from the manuscript about forward
% modeling coherent and incoherent neural sources to MEG responses.

% This figure shows the MEG forward model based on coherent and incoherent
% predictions coming from vertices located in V1 when using an absolute gain matrix.

% To runs this script, you need: 
% (1) Access to the SSMEG folder in the brainstorm data base
% (2) MEG_utils and Fieldtrip toolbox added to the paths. For example:
%     tbUse('ForwardModelSynchrony');
%        or to only add the MEG_utils toolbox:
%     addpath(genpath('~/matlab/git/toolboxes/meg_utils'))

%% 0. Set up paths and define parameters
figureDir       = fullfile(fmsRootPath, 'figures'); % Where to save images?
saveFigures     = true;     % Save figures in the figure folder?

% Path to brainstorm database
bs_db = '/Volumes/server/Projects/MEG/brainstorm_db/';

% Define project name, subject and data/anatomy folders
project_name = 'SSMEG';

% How many iterations of smoothing??
iterations = 'phase_0'; % pick between 'phase_1', 'phase_5', 'phase_10', 'phase_20', 'phase_nosmoothing', 'phase_1_normalized' where the number stands for the number of smoothing iterations

% Which subjects to average?
subject = {'wl_subj002','wl_subj004','wl_subj005','wl_subj006','wl_subj010','wl_subj011'};


%% Load subject predictions

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
    G_constrained = abs(G_constrained);
    
    % Load V1-3 template in downsampled brainstorm format (computed by interp_retinotopy.m)
    template = load(fullfile(anat_dir, 'areas_overlay.mat')); % [1xNsources] Every value between [-3,3] is inside V1-3, zeros refer to outside of visual cortex
    V1template = abs(template.sub_bs_areas)==1;
    
    V123PhaseScrambledTemplate = load(fullfile(anat_dir, sprintf('areas_overlay_%s.mat',iterations))); % [1xNsources] Every value between [-3,3] is inside V1-3, zeros refer to outside of visual cortex
    
    % Load V1-3 eccentricity in downsampled brainstorm format
    eccen    = load(fullfile(anat_dir, 'eccen_overlay.mat')); % [1xNsources] Every value from [1 3] is inside V1-3, zeros refer to outside of visual cortex
    polarang = load(fullfile(anat_dir, 'angle_overlay.mat'));
    polarang.sub_bs_angle_rad = pi/180*(90-polarang.sub_bs_angle);
    
    % Restrict V1 to 11 degrees eccentricity (stimulus diameter = 22)
    V1templateStimEccen = V1template.*(eccen.sub_bs_eccen<=11);
    V1PhaseScrambledTemplateStimEccen   = V123PhaseScrambledTemplate.sub_bs_areas.*V1templateStimEccen;
    
    %% 2. Compute forward solution (Method 1)
    % Compute the sensor weights, w, from V1-V3 using the contrained gain field (forward model)
%     w = G_constrained*V1template'; %  Nsensors x 1;
    
    phAmp2complex = @(r,th) r .* exp(1i*th);
    % Make a complex number with amplitudes and scrambled phases:
    V1complex = phAmp2complex(V1templateStimEccen,V1PhaseScrambledTemplateStimEccen);    
    
    wComplexV1(s,:) = G_constrained*V1complex';

    % Compute sensor weights for restricted by stimulus eccentricity
    w1_stimsize(s,:) = G_constrained*V1templateStimEccen'; %  Nsensors x 1;
    
    
    
end

% plot for subj 002
fH1 = figure(1); clf; set(fH1, 'Position', [1 1 1600 513], 'Name', 'Figure 3: V1 model predictions');
% rg = [-1 1]*10E-5;
rg = 0.5*[-1 1];

idx = isfinite(w1_stimsize(1,1:157));
absnorm = @(x, idx) abs(x) ./ norm(abs(x(idx)));


ch = subplot(221); 
megPlotMap(absnorm(w1_stimsize(1,1:157),idx),rg,[],bipolar,[],[],[],'isolines', 1);
set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12); title('Uniform phase S1')

ch = subplot(222); 
megPlotMap(absnorm(wComplexV1(1,1:157),idx),rg,[],bipolar,[],[],[],'isolines', 1)
set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12); title('Random phase S1')

%% Compute the mean across subjects
mnScrambledPhaseV1 = mean(wComplexV1,1);
mnUniformPhaseV1 = mean(w1_stimsize,1);

idx = isfinite(mnUniformPhaseV1(1:157));

ch = subplot(223); 
megPlotMap(absnorm(mnUniformPhaseV1(1:157),idx),rg,[],bipolar,[],[],[],'isolines', 1)
set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12); title('Uniform phase Average (S1-S6)')

ch = subplot(224); 
megPlotMap(absnorm(mnScrambledPhaseV1(1:157),idx),rg,[],bipolar,[],[],[],'isolines', 1)
set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12); title('Random phase Average (S1-S6)')

if saveFigures
    hgexport(gcf,fullfile(figureDir, ['Figure5_predictionV1_oneVSaverage_absGain.eps']));
%     figurewrite(fullfile(figureDir, ['predictionV1_oneVSaverage']), [],0,'.',1);
end
