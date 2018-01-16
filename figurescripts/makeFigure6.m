function makeFigure6()

% This function is under construction

% This function uses the forward model to create weights for each timepoint
% of the modeled MEG data. Weights are then put into sqd file, which needs 
% to be imported into the database via Brainstorm GUI. So far, this has
% only been done for wl_subj010.

% Once in Brainstorn, an inverse soluation can be made. We use this inverse
% solution to create a model of time varying sources. 



% Question:
% Should we use the initial headmodel or the new headmodel that was created
% based on the time varying weights?


%% 0. Set up paths and define parameters

% Brainstorm Database path
bs_db = '/Volumes/server/Projects/MEG/brainstorm_db/';

% Define project name, subject and data/anatomy folders
project_name = 'SSMEG';
subject = 'wl_subj010'; % pick 02, 04, 05, 06, 10, 11
iterations = 'phase_0'; % number stands for the number of smoothing iterations

d = dir(fullfile(bs_db, project_name, 'data', subject));
if strcmp(subject,'wl_subj002')
    data_dir = fullfile(bs_db, project_name, 'data', subject, d(6).name);
else
    data_dir = fullfile(bs_db, project_name, 'data', subject, d(5).name);
end

anat_dir = fullfile(bs_db, project_name, 'anat', subject);


%% 1. Create weights of forwardmodel

% Load headmodel from Brainstorm
headmodel = load(fullfile(data_dir, 'headmodel_surf_os_meg.mat'));

% Get Gain matrix
G = headmodel.Gain; % [Nsensors x 3*Nvertices]

% Contrained gain matrix
G_constrained = bst_gain_orient(G, headmodel.GridOrient); % [Nsensors x Nsources], equivalent to size BS pial cortex [1x15002]

% Load V1-3 template with unitized phases in downsampled brainstorm format (computed by interp_retinotopy.m)
areas = load(fullfile(anat_dir, 'areas_overlay.mat')); % [1xNsources] Every value between [-3,3] is inside V1-3, zeros refer to outside of visual cortex
eccen    = load(fullfile(anat_dir, 'eccen_overlay.mat')); % [1xNsources] Every value from [1 3] is inside V1-3, zeros refer to outside of visual cortex
polarang = load(fullfile(anat_dir, 'angle_overlay.mat'));

% Get only vertices in V1
template.V1     = abs(areas.sub_bs_areas)==1;

% Get radians instead of angles
polarang.sub_bs_angle_rad = pi/180*(90-polarang.sub_bs_angle);

% Limit to 11 degrees eccentricity
template.V1StimEccenAmplitudes = template.V1.*(eccen.sub_bs_eccen<=11);

%% 2. Get weights for every time point

% Create a 12 Hz sine wave
dt = .01;        % s
t  = 0:dt:.24; % s
f  = 12;          % Hz

sin12 = sin(2*pi*f*t);
% sin12 = reshape(sin12, [1 size(sin12)]);

phAmp2complex = @(r,th) r .* exp(1i*th);

% Create sources with the same amplitude (value of 1), but either incoherent or coherent phase
for ii = 1:length(t)
    
    % for coherent phase
    thisPhase = template.V1StimEccenAmplitudes .* sin12(ii);
    template.V1StimEccenPhaseC(:,:,ii) = phAmp2complex(template.V1StimEccenAmplitudes,thisPhase);
    
    % for incoherent phase
    thisPhase = template.V1StimEccenAmplitudes .* (2*pi*rand(size(template.V1StimEccenAmplitudes)));
    template.V1StimEccenPhaseI(:,:,ii) = phAmp2complex(template.V1StimEccenAmplitudes,thisPhase);
end
    
% Get weights 
for ii = 1:size(template.V1StimEccenPhaseC,3)    
    w.V1c(:,ii) = squeeze(G_constrained*template.V1StimEccenPhaseC(:,:,ii)'); %  Nsensors x time point;
    w.V1i(:,ii) = squeeze(G_constrained*template.V1StimEccenPhaseI(:,:,ii)'); %  Nsensors x time points;
end

% Visualize

% Get colorbar limits
clims = [-1,1] .* abs(max(max(w.V1c)));

% Visualize weights
figure; 
for ii = 1:length(t)
    clf;
    megPlotMap(abs(w.V1c(1:157,ii)),clims,[],bipolar,sprintf('Forward model for timepoint %d',ii));
    pause(.1);
end

figure; 
for ii = 1:length(t)
    clf;
    megPlotMap(abs(w.V1i(1:157,ii)),clims,[],bipolar,sprintf('Forward model for timepoint %d',ii));
    pause(.1);
end

figure; megPlotMap(mean(abs(w.V1c(1:157,:)),2),clims,[],bipolar,'Mean weights across timepoints: Coherent')
figure; megPlotMap(mean(abs(w.V1i(1:157,:)),2),clims,[],bipolar,'Mean weights across timepoints: Incoherent')


%% 3. Create SQD file (Only necessary once, after that it should be saved in Brainstorm DB) 

% [~, meg_files] = meg_load_sqd_data('/Volumes/server/Projects/MEG/SSMEG/09_SSMEG_06_27_2014_wl_subj010/raw/','V1ForwardStimEccen_12HzCycle');
% 
% newFile = '~/Desktop/testV1Forward12Hz.sqd';
% newFile = '~/Desktop/testV1ForwardRandom.sqd';
% 
% % new file needs to be imported into Brainstorm session
% sqdwrite(fullfile(meg_files.folder,meg_files.name),newFile,w.V1c');
% sqdwrite(fullfile(meg_files.folder,meg_files.name),newFile,w.V1i');


%% 4. Create inverse model

% Load Brainstorm downsampled pial surface
bs_pial_low = load(fullfile(anat_dir, 'tess_cortex_pial_low.mat'));

% Load Brainstorm inverse model
inverseC  = load(fullfile(bs_db, project_name, 'data', subject, 'testV1Forward12Hz/results_MN_MEG_KERNEL_180116_1832.mat'));
inverseI  = load(fullfile(bs_db, project_name, 'data', subject, 'testV1ForwardRandom/results_MN_MEG_KERNEL_180116_1828.mat'));

% Create source response for each timepoint
for ii = 1:size(w.V1c,2)
    s.V1c(:,ii) = inverseC.ImagingKernel*w.V1c(1:157,ii);
    s.V1i(:,ii) = inverseI.ImagingKernel*w.V1i(1:157,ii);
end


s = struct2cell(s);

for source = 1:size(s,1)
    
    thisSource = s{source};

    % Visualize: set up mesh
    figure; tH = trimesh(bs_pial_low.Faces,bs_pial_low.Vertices(:,1),bs_pial_low.Vertices(:,2),bs_pial_low.Vertices(:,3));
    
    % Visualize: set curvature colors
    colors = zeros(size(bs_pial_low.Vertices,1),1);
    colors(bs_pial_low.Curvature<0) = -1.5;
    colors(bs_pial_low.Curvature>=0) = -.5;
    
    % Define colors
    cmap = [gray(128); jet(128)];
    axis equal; hold on
    
    % Plot it
    for ii = 1:size(thisSource,2)
        
        colors = abs(thisSource(:,ii));
        
        % set source pediction as colors
        set(tH, 'LineStyle', 'none', 'FaceColor', 'interp', 'FaceVertexCData',colors);
        drawnow;
        colormap(cmap); set(gca, 'CLim',2*pi*[-.01 .01]); colorbar;
        
        pos = [-.1 0 .1];
        light('Position',pos,'Style','local')
        lighting gouraud
        material shiny; %dull
        title(ii);
        pause(.1);
    end
    
    % Plot mean timeseries of V1 sources
    timeseriesV1 = thisSource(logical(template.V1StimEccenAmplitudes'),:);
    figure; plot(t,mean(abs(timeseriesV1),1))
    
    % Plot mean timeseries of non visual sources
    nonVisualVertices = abs(areas.sub_bs_areas)==0;
    figure; plot(t,mean(thisSource(nonVisualVertices,:),1));

end

%% NOT READY YET: Save to FS space

% BRAINSTORM GUI HAS TO BE OPEN FOR THIS STEP

% if ~exist('MRIwrite')
%     addpath(genpath('/Applications/freesurfer/matlab'));
%     addpath(genpath('/Applications/freesurfer/fsfast/toolbox'));
% % addpath(genpath('~/matlab/git/toolboxes/brainstorm3'));
% end
% 
% 
% % Where are the FS subjects?
% fs_dir = '/Volumes/server/Freesurfer_subjects/';
% 
% % Set Brainstorm mesh back to FS mesh
% [fs_lh_overlay fs_rh_overlay] = tess_bst2fs('wl_subj010', fullfile(fs_dir,'wl_subj010'), s.V1);
% 
% 
% MRIwrite(struct('vol', fs_lh_overlay), fullfile(fs_dir,'wl_subj010','surf','lh.inverse_V1_coherent.mgz'));
% MRIwrite(struct('vol', fs_rh_overlay), fullfile(fs_dir,'wl_subj010','surf','rh.inverse_V1_coherent.mgz'));

return