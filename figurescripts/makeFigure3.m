function makeFigure3()

% This is a function to make Figure 3 from the manuscript about forward
% modeling coherent and incoherent neural sources to MEG responses.

% This figure shows the MEG forward model based on coherent and incoherent
% predictions coming from vertices located in V1.

% To runs this script, you need:     
% (1) Access to the SSMEG folder in the brainstorm data base
% (2) MEG_utils and Fieldtrip toolbox added to the paths. For example:
%     tbUse('ForwardModelSynchrony');
%        or to only add the MEG_utils toolbox:
%     addpath(genpath('~/matlab/git/toolboxes/meg_utils'))


% Set up paths
figureDir       = fullfile(fmsRootPath, 'figures'); % Where to save images?
saveFigures     = true;     % Save figures in the figure folder?


%% 0. Set up paths and define parameters
% Path to brainstorm database
bs_db = '/Volumes/server/Projects/MEG/brainstorm_db/';

% Define project name, subject and data/anatomy folders
project_name = 'SSMEG';

% Which subjects to average?
subject = {'wl_subj002','wl_subj004','wl_subj005','wl_subj006','wl_subj010','wl_subj011'};

% What's the plotting range
climsOne = [-1 1]*1E-4;
climsAve = [-.5 .5]*1E-4;

% Number of iterations for the random coherence prediction of the forward
% model
iter = 1000;

for s = 1:length(subject)
    
    d = dir(fullfile(bs_db, project_name, 'data', subject{s}, 'R*'));
    data_dir = fullfile(d(1).folder, d(1).name);    
    anat_dir = fullfile(bs_db, project_name, 'anat', subject{s});
    
    %% 1. Load relevant matrices
    % Load Gain matrix created by brainstorm
    headmodel = load(fullfile(data_dir, 'headmodel_surf_os_meg.mat'));
    G = headmodel.Gain(1:157,:); % Gain matrix, [Nsensors x 3*Nvertices]
    G_constrained = bst_gain_orient(G, headmodel.GridOrient); % Contrained gain matrix [Nsensors x Nsources], equivalent to size pial cortex [1x15002]
    
    % Load V1-3 template in downsampled brainstorm format (computed by interp_retinotopy.m)
    areas = load(fullfile(anat_dir, 'areas_overlay.mat')); % [1xNsources] Every value between [-3,3] is inside V1-3, zeros refer to outside of visual cortex
    template.V1 = abs(areas.sub_bs_areas)==1;
    
%     V123PhaseScrambledTemplate = load(fullfile(anat_dir, sprintf('areas_overlay_%s.mat',iterations))); % [1xNsources] Every value between [-3,3] is inside V1-3, zeros refer to outside of visual cortex
    
    % Load V1-3 eccentricity in downsampled brainstorm format
    eccen    = load(fullfile(anat_dir, 'eccen_overlay.mat')); % [1xNsources] Every value from [1 3] is inside V1-3, zeros refer to outside of visual cortex
    polarang = load(fullfile(anat_dir, 'angle_overlay.mat'));
    polarang.sub_bs_angle_rad = pi/180*(90-polarang.sub_bs_angle);
    
    % Restrict V123 to 11 degrees eccentricity (stimulus diameter = 22)
    template.V1StimEccenAmplitudes = template.V1 .*(eccen.sub_bs_eccen<=11);
    
    template.V1StimEccenPhaseCoherent     = template.V1StimEccenAmplitudes * 0;    
    template.V1StimEccenPhaseIncoherent   = template.V1StimEccenAmplitudes .* (rand(iter,size(template.V1StimEccenAmplitudes,2)) * 2*pi);    

    %% 2. Compute forward solution (Method 1)
    % Compute the sensor weights, w, from V1-V3 using the contrained gain field (forward model)
%     w = G_constrained*V1template'; %  Nsensors x 1;
    
    phAmp2complex = @(r,th) r .* exp(1i*th);
    
    % Make a complex number with amplitudes and phases:
    template.V1coherent = phAmp2complex(template.V1StimEccenAmplitudes,template.V1StimEccenPhaseCoherent);
     
    template.V1incoherent = phAmp2complex(repmat(template.V1StimEccenAmplitudes,[iter,1]),template.V1StimEccenPhaseIncoherent);

    % Compute prediction of forward model for template restricted by stimulus eccentricity
    w.V1c(s,:) = G_constrained*template.V1coherent'; %  Nsensors x 1;
    w.V1i(s,:) = mean(abs(G_constrained*template.V1incoherent'),2); %  Nsensors x 1;
    
end

%% 3. Visualize predictions from forward model

% plot for subj 002
fH1 = figure(1); clf; set(fH1, 'Position', [1 1 1600 513], 'Name', 'Figure 3: V1 model predictions');

ch = subplot(221); 
megPlotMap(abs(w.V1c(1,:)),climsOne,[],bipolar,[],[],[],'isolines', 0.75*max(climsOne) * [1 1]);
set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12); title('Coherent phase S1')

ch = subplot(222); 
megPlotMap(abs(w.V1i(1,:)),climsOne,[],bipolar,[],[],[],'isolines',  0.75*max(climsOne) * [1 1])
set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12); title('Incoherent phase S1')

%% Compute the mean across subjects

w.V1c_mn = mean(w.V1c,1);
w.V1i_mn = mean(w.V1i,1);

ch = subplot(223); 
megPlotMap(abs(w.V1c_mn),climsAve,[],bipolar,[],[],[],'isolines', 0.75*max(climsAve) * [1 1])
set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12); title('Coherent phase Average (S1-S6)')

ch = subplot(224); 
megPlotMap(abs(w.V1i_mn),climsAve,[],bipolar,[],[],[],'isolines',  0.75*max(climsAve) * [1 1])
set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12); title('Incoherent phase Average (S1-S6)')

if saveFigures
    hgexport(gcf,fullfile(figureDir, 'Figure3_predictionV1_oneVSaverage.eps'));
%     figurewrite(fullfile(figureDir, ['predictionV1_oneVSaverage']), [],0,'.',1);
end
