function makeFigure3()

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
%     V123complex = phAmp2complex(V123templateStimEccen,V123PhaseScrambledTemplateStimEccen);
    V1complex = phAmp2complex(V1templateStimEccen,V1PhaseScrambledTemplateStimEccen);
%     V2complex = phAmp2complex(V2templateStimEccen,V2PhaseScrambledTemplateStimEccen);
%     V3complex = phAmp2complex(V3templateStimEccen,V3PhaseScrambledTemplateStimEccen);
    
    
%     wComplexV123(s,:) = G_constrained*V123complex';
    wComplexV1(s,:) = G_constrained*V1complex';
%     wComplexV2(s,:) = G_constrained*V2complex';
%     wComplexV3(s,:) = G_constrained*V3complex';
    
    % Compute sensor weights for restricted by stimulus eccentricity
%     w_stimsize(s,:)  = G_constrained*V123templateStimEccen'; %  Nsensors x 1;
    w1_stimsize(s,:) = G_constrained*V1templateStimEccen'; %  Nsensors x 1;
%     w2_stimsize(s,:) = G_constrained*V2templateStimEccen'; %  Nsensors x 1;
%     w3_stimsize(s,:) = G_constrained*V3templateStimEccen'; %  Nsensors x 1;
    
    
    
end

% plot for subj 002
fH1 = figure(1); clf; set(fH1, 'Position', [1 1 1600 513], 'Name', 'Figure 3: V1 model predictions');
rg = [-1 1]*10E-5;

ch = subplot(221); 
megPlotMap(abs(w1_stimsize(1,1:157)),rg,[],bipolar,[],[],[],'isolines', 1);
set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12); title('Uniform phase S2')

ch = subplot(222); 
megPlotMap(abs(wComplexV1(1,1:157)),0.5*rg,[],bipolar,[],[],[],'isolines', 1)
set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12); title('Random phase S2')

%% Compute the mean across subjects
% mnScrambledPhaseV123 = mean(wComplexV123,1);
mnScrambledPhaseV1 = mean(wComplexV1,1);
% mnScrambledPhaseV2 = mean(wComplexV2,1);
% mnScrambledPhaseV3 = mean(wComplexV3,1);

% mnUniformPhaseV123 = mean(w_stimsize,1);
mnUniformPhaseV1 = mean(w1_stimsize,1);
% mnUniformPhaseV2 = mean(w2_stimsize,1);
% mnUniformPhaseV3 = mean(w3_stimsize,1);


ch = subplot(223); 
megPlotMap(abs(mnUniformPhaseV1(1:157)),rg,[],bipolar,[],[],[],'isolines', 1)
set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12); title('Uniform phase Average (n=6)')

ch = subplot(224); 
megPlotMap(abs(mnScrambledPhaseV1(1:157)),0.5*rg,[],bipolar,[],[],[],'isolines', 1)
set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12); title('Random phase Average (n=6)')

if saveFigures
    hgexport(gcf,fullfile(figureDir, ['Figure3_predictionV1_oneVSaverage.eps']));
%     figurewrite(fullfile(figureDir, ['predictionV1_oneVSaverage']), [],0,'.',1);
end
