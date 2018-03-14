% Create and save meshes
if ~exist('fast_fileexists')
    addpath(genpath('/Applications/freesurfer/matlab'));
    addpath(genpath('/Applications/freesurfer/fsfast/toolbox'));
end


pth = '/Volumes/server/Freesurfer_subjects/wl_subj004/surf';
hemi = 'l';
stimEccen = 11; % stimulus had 11 degrees of visual angle


% Surface mesh
surf = fullfile(pth, sprintf('%sh.pial', hemi));
[vertices, faces] = freesurfer_read_surf(surf);

% Curvature (for coloring)
fname = fullfile(pth, sprintf('%sh.curv', hemi));
curv = read_curv(fname);
mshcolor.curv = zeros(size(curv));
mshcolor.curv(curv>0) = 1;

% Surface maps
%       angle, eccentrity, and area from Benson et al
retEccen = fullfile(pth, sprintf('%sh.template_eccen.mgz', hemi));
retArea  = fullfile(pth, sprintf('%sh.template_areas.mgz', hemi));
% retV1label = fullfile(pth, sprintf('%sh.template_V1.label', hemi));
ni       = MRIread(retEccen); mshcolor.eccen = ni.vol(:)<stimEccen;
ni       = MRIread(retArea);  mshcolor.areas = abs(ni.vol(:));
% ni       = MRIread(retV1label);  mshcolor.v1label = abs(ni.vol(:));


% Get V1 and restrict to stimulus eccentricity
mshcolor.v1 = mshcolor.areas==1;
mshcolor.v1StimEccen = mshcolor.v1.*mshcolor.eccen;
idx      = mshcolor.v1StimEccen > 0;

cmap = parula(256); clim = [0 1];

lightColor = [0.5000 0.5000 0.4500];
cmap(1,:) =  2 *lightColor; % color for gyri
cmap(2,:) =  1 *lightColor; % color for sulci



%% Plot

fH = figure; pos = get(gcf, 'Position'); pos([3 4]) = [1000 800];
set(fH, 'Color', 'w', 'Position', pos); clf;

% Mesh
c = mshcolor.curv;
c(idx) = round(mshcolor.v1StimEccen(idx)/clim(2)*60);

mx = vertices(:,1);%+surf_offsets(1);
my = vertices(:,2);%+surf_offsets(2);
mz = vertices(:,3);%+surf_offsets(3);

tH = trimesh(faces, mx, my, mz, c); axis equal; hold on
set(tH, 'LineStyle', 'none', 'FaceColor', 'interp', 'FaceVertexCData',c)
colormap(cmap); set(gca, 'CLim', [0 255]);

axis off;

% Viewpoint
if strcmpi(hemi, 'r'), set(gca, 'View', [-70 -10]);
else, set(gca, 'View', [90.2000, -6.2000]); end

% Lighting
if strcmpi(hemi, 'r'), pos = [-1 1 1]; else, pos = [1 1 1]; end
light('Position',100*pos,'Style','local')
lighting gouraud
material dull

