% Create and save meshes
if ~exist('fast_fileexists')
    addpath(genpath('/Applications/freesurfer/matlab'));
    addpath(genpath('/Applications/freesurfer/fsfast/toolbox'));
end


pth = '/Volumes/server/Freesurfer_subjects/wlsubj005/surf';
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
retAngle = fullfile(pth, sprintf('%sh.template_angle.mgz', hemi));
retEccen = fullfile(pth, sprintf('%sh.template_eccen.mgz', hemi));
retArea  = fullfile(pth, sprintf('%sh.template_areas.mgz', hemi));

% Angle:
ni             = MRIread(retAngle);
mshcolor.angle = ni.vol(:);


% Eccen:
ni                 = MRIread(retEccen);
mshcolor.stimEccen = ni.vol(:)<stimEccen;
mshcolor.eccen     = ni.vol(:);

% Area:
ni             = MRIread(retArea); 
mshcolor.areas = abs(ni.vol(:));


% Get V1 and restrict to stimulus eccentricity
mshcolor.v1 = mshcolor.areas==1;
mshcolor.v123 = mshcolor.areas>0;

mshcolor.v1StimEccen = mshcolor.v1.*mshcolor.stimEccen;
mshcolor.v123StimEccen = mshcolor.v123.*mshcolor.stimEccen;

idx      = mshcolor.v123StimEccen > 0;

cmap = parula(256); clim = [0 1];

lightColor = [0.5000 0.5000 0.4500];
cmap(1,:) =  2 *lightColor; % color for gyri
cmap(2,:) =  1 *lightColor; % color for sulci

grn = cmap(150,:);

%% Plot V1-V3

fH = figure; pos = get(gcf, 'Position'); pos([3 4]) = [1000 800];
set(fH, 'Color', 'w', 'Position', pos); clf;

% Mesh
c = mshcolor.curv;
c(idx) = mshcolor.v123StimEccen(idx).*150;

mx = vertices(:,1);%+surf_offsets(1);
my = vertices(:,2);%+surf_offsets(2);
mz = vertices(:,3);%+surf_offsets(3);

tH = trimesh(faces, mx, my, mz, c); axis equal; hold on
set(tH, 'LineStyle', 'none', 'FaceColor', 'interp', 'FaceVertexCData',c)
colormap(cmap); set(gca, 'CLim', [0 255]);

axis off;

% Viewpoint
if strcmpi(hemi, 'r'), set(gca, 'View', [-90 -4]);
else, set(gca, 'View', [90.2000, -6.2000]); end

% Lighting
if strcmpi(hemi, 'r'), pos = [-1 1 1]; else, pos = [1 1 1]; end
light('Position',100*pos,'Style','local')
lighting gouraud
material dull


%% Plot Eccen

fH = figure; pos = get(gcf, 'Position'); pos([3 4]) = [1000 800];
set(fH, 'Color', 'w', 'Position', pos); clf;

% Mesh
idx = mshcolor.eccen > 0;
c = mshcolor.curv;
c(c==1)= 1/256;
c(idx) = 2*(1/256) + mshcolor.eccen(idx)./90;

% Get xyz vertices
mx = vertices(:,1);%+surf_offsets(1);
my = vertices(:,2);%+surf_offsets(2);
mz = vertices(:,3);%+surf_offsets(3);

% get colors correctly:
lightColor = [0.5000 0.5000 0.4500];
sulci = 1 *lightColor;
gyri  = 2 *lightColor;
cmap = [gyri; sulci; hsv(256)]; clim = [0 1];

% cmap(1,:) =  2 *lightColor; % color for gyri
% cmap(2,:) =  1 *lightColor; % color for sulci


tH = trimesh(faces, mx, my, mz, c); axis equal; hold on
set(tH, 'LineStyle', 'none', 'FaceColor', 'interp', 'FaceVertexCData',c)
colormap(cmap); set(gca, 'CLim', clim);

axis off;

% Viewpoint
if strcmpi(hemi, 'r'), set(gca, 'View', [-90 -4]);
else, set(gca, 'View', [90.2000, -6.2000]); end

% Lighting
if strcmpi(hemi, 'r'), pos = [-1 1 1]; else, pos = [1 1 1]; end
light('Position',100*pos,'Style','local')
lighting gouraud
material dull

%% Plot Angle

fH = figure; pos = get(gcf, 'Position'); pos([3 4]) = [1000 800];
set(fH, 'Color', 'w', 'Position', pos); clf;

% Mesh
idx = mshcolor.angle ~= 0;
c = mshcolor.curv;
c(c==1)= 1/34;
c(idx) = (2/34) + mshcolor.angle(idx)./180;

% Get xyz vertices
mx = vertices(:,1);%+surf_offsets(1);
my = vertices(:,2);%+surf_offsets(2);
mz = vertices(:,3);%+surf_offsets(3);

% get colors correctly:
lightColor = [0.5000 0.5000 0.4500];
sulci = 1 *lightColor;
gyri  = 2 *lightColor;
cmapLR = cmapangLR;
if hemi == 'l'
    cmapLR = cmapLR([(1:16),((64-13):64)],:);
else
     cmapLR = cmapLR((17:48),:);
end
cmap = [gyri; sulci; cmapLR]; clim = [0 1];

% cmap(1,:) =  2 *lightColor; % color for gyri
% cmap(2,:) =  1 *lightColor; % color for sulci


tH = trimesh(faces, mx, my, mz, c); axis equal; hold on
set(tH, 'LineStyle', 'none', 'FaceColor', 'interp', 'FaceVertexCData',c)
colormap(cmap); set(gca, 'CLim', clim);

axis off;

% Viewpoint
if strcmpi(hemi, 'r'), set(gca, 'View', [-90 -4]);
else, set(gca, 'View', [90.2000, -6.2000]); end

% Lighting
if strcmpi(hemi, 'r'), pos = [-1 1 1]; else, pos = [1 1 1]; end
light('Position',100*pos,'Style','local')
lighting gouraud
material dull

%% Plot 11 degree aperture

fH = figure; pos = get(gcf, 'Position'); pos([3 4]) = [1000 800];
set(fH, 'Color', 'w', 'Position', pos); clf;

% Mesh
idx1 = ((0 < mshcolor.eccen)&(mshcolor.eccen<= 11));
idx2 = ((11 < mshcolor.eccen)&(mshcolor.eccen<= 90));

c = mshcolor.curv;
c(c==1) = 1/4;
c(idx1) = 2/3;
c(idx2) = 1;

% Get xyz vertices
mx = vertices(:,1);%+surf_offsets(1);
my = vertices(:,2);%+surf_offsets(2);
mz = vertices(:,3);%+surf_offsets(3);

% get colors correctly:
lightColor = [0.5000 0.5000 0.4500];
sulci = 1 *lightColor;
gyri  = 2 *lightColor;
whte  = [1 1 1];
blk   = [0.25 0.25 0.25];
cmap = [gyri; sulci;  grn;  blk]; clim = [0 1];

% cmap(1,:) =  2 *lightColor; % color for gyri
% cmap(2,:) =  1 *lightColor; % color for sulci


tH = trimesh(faces, mx, my, mz, c); axis equal; hold on
set(tH, 'LineStyle', 'none', 'FaceColor', 'interp', 'FaceVertexCData',c)
colormap(cmap); set(gca, 'CLim', clim);

axis off;

% Viewpoint
if strcmpi(hemi, 'r'), set(gca, 'View', [-90 -4]);
else, set(gca, 'View', [90.2000, -6.2000]); end

% Lighting
if strcmpi(hemi, 'r'), pos = [-1 1 1]; else, pos = [1 1 1]; end
light('Position',100*pos,'Style','local')
lighting gouraud
material dull

%% Make colormap 
% figure(1); clf;
% R = 80:-0.3125:0.1; % radii
% N = numel(R) ; % 
% Cmap = fliplr(hsv(N)')' ; % one of the many available maps, see the help
% % now loop over the radii and colormaps using indexing
% for k = 1:N
%     pos = R(k) * [-1 -1 2 2] ;    % position
%     clr = Cmap(k,:) ;             % RGB triplet
%     rectangle('Position', pos, 'Curvature', [1 1], ... 
%               'EdgeColor',clr, 'LineWidth', 4)
% end
% axis square

% figure(2); clf;
% r = linspace(0,1,10);
% theta = linspace(0.5*pi, -1.5*pi, 100);
% [rg, thg] = meshgrid(r,theta);
% [x,y] = pol2cart(thg,rg);
% pcolor(x,y,thg);
% colormap(hsv);
% shading flat;
% axis equal; axis off

figure; drawcolorbarcircular(cmapangLR,1)