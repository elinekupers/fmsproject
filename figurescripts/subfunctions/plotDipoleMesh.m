% Brainstorm Database path
bsDB  = '/Volumes/server/Projects/MEG/brainstorm_db/';

% Freesurfer subject path
fsDir = '/Volumes/server/Freesurfer_subjects/';

% Define project name, subject and data/anatomy folders
project_name    = 'SSMEG';
subject         = 'wl_subj010'; % pick 02, 04, 05, 06, 10, 11

% Get directories
d               = dir(fullfile(bsDB, project_name, 'data', subject, 'R*'));
dataDir         = fullfile(d(1).folder, d(1).name);
anatDir         = fullfile(bsDB, project_name, 'anat', subject);

bs_pial_low = load(fullfile(anatDir, 'tess_cortex_pial_low.mat'));

% Get V1 template limited to 11 degrees eccentricity
template = getTemplate(anatDir, 'V1', 11);

colors = abs(template.V1StimEccen);
    
% Define colorbar colors
cmap = [gray(128); jet(128)];

% Plot mean amplitude across epochs
figure; set(gcf, 'Color', 'w', 'Position', [163 483 891 554])

% Set up mesh
tH = trimesh(bs_pial_low.Faces,bs_pial_low.Vertices(:,1),bs_pial_low.Vertices(:,2),bs_pial_low.Vertices(:,3));
axis equal; hold on


quiver3(bs_pial_low.Vertices(:,1),bs_pial_low.Vertices(:,2),bs_pial_low.Vertices(:,3), ...
bs_pial_low.VertNormals(:,1),bs_pial_low.VertNormals(:,2),bs_pial_low.VertNormals(:,3))


% set source pediction as colors
set(tH, 'LineStyle', 'none', 'FaceColor', 'interp', 'FaceVertexCData',double(colors));

colormap(cmap); colorbar; set(gca, 'CLim', clims);

pos = [-.1 0 .1];
light('Position',pos,'Style','local')
lighting gouraud
material shiny; %dull

drawnow;
pause(.1);