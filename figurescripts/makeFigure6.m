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
template.V1StimEccen = template.V1.*(eccen.sub_bs_eccen<=11);

%% 2. Get weights for every time point

% Create a 12 Hz sine wave
dt = .001;        % s
t  = (0:dt:.999); % s
f  = 12;          % Hz

sin12 = sin(2*pi*f*t);
sin12 = reshape(sin12, [1 size(sin12)]);

template.V1StimEccenFullCycle = repmat(template.V1StimEccen, [1 1 length(sin12)]) .* sin12;

for t = 1:size(template.V1StimEccenFullCycle,3)    
    w.V1(:,:,t) = G_constrained*template.V1StimEccenFullCycle(:,:,t)'; %  Nsensors x 1;    
end

% Squeeze (probably can be omitted when cleaning up code and creating w.V1 as a matrix)
w.V1 = squeeze(w.V1);

% Get colorbar limits
clims = [min(min(w.V1)), max(max(w.V1))];

% Visualize weights
figure; clf;
for ii = 1:size(w.V1,2)/100 % only show first 
    megPlotMap(w.V1(1:157,ii),clims,[],bipolar,sprintf('Forward model for timepoint %d',ii));
    pause(.1);
end

figure; megPlotMap(mean(abs(w.V1(1:157,:)),2),clims,[],bipolar,'Mean weights across timepoints')

%% 3. Create SQD file (Only necessary once, after that it should be saved in Brainstorm DB) 

% [~, meg_files] = meg_load_sqd_data('/Volumes/server/Projects/MEG/SSMEG/09_SSMEG_06_27_2014_wl_subj010/raw/','V1ForwardStimEccen_12HzCycle');
% 
% newFile = '~/Desktop/testV1Forward12Hz.sqd';
% 
% % new file needs to be imported into Brainstorm session
% sqdwrite(fullfile(meg_files.folder,meg_files.name),newFile,squeeze(w.V1)');


%% 4. Create inverse model

% Load Brainstorm downsampled pial surface
bs_pial_low = load(fullfile(anat_dir, 'tess_cortex_pial_low.mat'));

% Load Brainstorm inverse model
inverse  = load(fullfile(bs_db, project_name, 'data', subject, 'testV1Forward12Hz/results_MN_MEG_KERNEL_180116_1206.mat'));

% Create source response for each timepoint
for t = 1:size(template.V1StimEccenFullCycle,3)
    s.V1(:,t) = inverse.ImagingKernel*w.V1(1:157,t);
end

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
for ii = 1:size(s.V1,2)
    colors = s.V1(:,ii);
    
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
timeseriesV1 = s.V1(logical(template.V1StimEccen'),:);
figure; plot(mean(timeseriesV1))

% Plot mean timeseries of non visual sources
nonVisualVertices = abs(areas.sub_bs_areas)==0;
figure; plot(mean(s.V1(nonVisualVertices,:),2));

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