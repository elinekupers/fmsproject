%% s_headmodelGeometryDemo.m

% Load subject data
subject         = 'wlsubj070';
bsDB            = '/Volumes/server/Projects/MEG/brainstorm_db/';
projectName     = 'SSMEG';
keep_sensors    = logical([ones(157,1); zeros(192-157,1)]); % Note: Figure out a more generic way to define keep_sensors

% Brainstorm mesh
d = dir(fullfile(bsDB, projectName, 'data', subject, 'R*'));
bsData = fullfile(d(1).folder, d(1).name);
bsAnat = fullfile(bsDB, projectName, 'anat', subject);

% fsAnat = fullfile(sprintf('/Volumes/server/Freesurfer_subjects/%s/surf/',subject));

%% 1. Load Gain matrix

G_constrained = getGainMatrix(bsData, keep_sensors);

%% 2. Load retinotopy templates (downsampled)

% Load V1-3 template with unitized phases in downsampled brainstorm format (computed by interp_retinotopy.m)
areas    = load(fullfile(bsAnat, 'areas_overlay.mat')); % [1xNsources] Every value between [-3,3] is inside V1-3, zeros refer to outside of visual cortex. CHECK: Positive values represent LH (?) Negative values RH (?)
eccen    = load(fullfile(bsAnat, 'eccen_overlay.mat')); % [1xNsources] Nonzero value represents vertex preferred eccentricity in degrees, zeros refer to outside of visual cortex
polarang = load(fullfile(bsAnat, 'angle_overlay.mat')); % [1xNsources] Values represents vertex preferred polar angle in degrees (0 is the upper vertical meridian), zeros refer to outside of visual cortex

% Get V1 template
template.V1       = abs(areas.sub_bs_areas)==1;

% For reference, get polar angle of each vertex in degrees (not used yet)
polarang.sub_bs_angle_rad = pi/180*(90-polarang.sub_bs_angle);

%% Get BS pial surface
bs_pial_low = load(fullfile(bsAnat, 'tess_cortex_pial_low.mat'));

%% Plot one vertex and its forward model projection
eccen = 4;
pa    = [88,88.2];

[~,idx] = find((polarang.sub_bs_angle>pa(1))&(polarang.sub_bs_angle<pa(2)));
thisVertex1  = idx(1);

pa    = [91.8,92];
[~,idx] = find((polarang.sub_bs_angle>pa(1))&(polarang.sub_bs_angle<pa(2)));
thisVertex2  = idx(1);

%% Plot opposite vertices predictions

figure(1); clf; 
subplot(211); megPlotMap(G_constrained(:,thisVertex1),1E-5.*[-0.5, 0.5],[],'bipolar', 'Prediction: Upper visual field vertex (PA 88 deg)');

subplot(212); megPlotMap(G_constrained(:,thisVertex2),1E-5.*[-1, 1],[],'bipolar','Prediction: Lower visual field vertex (PA 92 deg)');

% Define colorbar colors
cmap = [gray(128); jet(128)];
clims = [-2 2];

%% Plot opposite vertices in brainstorm mesh

figure(2); clf;
tH = trimesh(bs_pial_low.Faces,bs_pial_low.Vertices(:,1),bs_pial_low.Vertices(:,2),bs_pial_low.Vertices(:,3));
axis equal; hold on

curv = bs_pial_low.Curvature - max(bs_pial_low.Curvature);
curv(thisVertex1) = 1;
curv(thisVertex2) = 2;

alphaVertices = ones(size(bs_pial_low.Vertices,1),1).*0.1;
alphaVertices(thisVertex1) = .1;
alphaVertices(thisVertex2) = .1;

% set source pediction as colors
set(tH, 'LineStyle', 'none', 'FaceColor', 'flat', 'FaceVertexCData',double(curv), 'FaceAlpha', 'flat', 'FaceVertexAlphaData', alphaVertices);

colormap(cmap); colorbar; set(gca, 'CLim',clims);

pos = [-.1 0 .1];
light('Position',pos,'Style','local')
lighting gouraud
material shiny; %dull


%% Plot surface normals on brainstorm mesh

figure(3); clf; set(gcf, 'Color', 'w', 'Position', [163 483 891 554])

% Set up mesh
tH = trimesh(bs_pial_low.Faces,bs_pial_low.Vertices(:,1),bs_pial_low.Vertices(:,2),bs_pial_low.Vertices(:,3));
axis equal; hold on

quiver3(bs_pial_low.Vertices(:,1),bs_pial_low.Vertices(:,2),bs_pial_low.Vertices(:,3), ...
bs_pial_low.VertNormals(:,1),bs_pial_low.VertNormals(:,2),bs_pial_low.VertNormals(:,3), 2)

colormap(cmap); colorbar; set(gca, 'CLim',[0 1]);

