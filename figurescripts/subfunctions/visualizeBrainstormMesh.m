function visualizeBrainstormMesh(anatDir, data, thresh, cmap, clims, meshType, ttl)

%% Check inputs

if ~exist('meshType','var') || isempty(meshType)
    bs_pial_low = load(fullfile(anatDir, 'tess_cortex_pial_low.mat'));
else
    bs_pial_low = load(fullfile(anatDir, 'tess_cortex_pial_low.mat'));
    if strcmp(meshType,'smooth')
        smoothed = load(fullfile(anatDir, 'tess_cortex_pial_low_fig.mat'));
        bs_pial_low.Vertices    = smoothed.Vertices;
        bs_pial_low.Faces       = smoothed.Faces;
        bs_pial_low.Comment     = smoothed.Comment;
    end
end


if ~exist('cmap','var') || isempty(cmap)
    cmap = hsv(256);
    
end


if ~exist('thresh','var') || isempty(thresh)
    thresh = 0;
end

if ~exist('clims','var') || isempty(clims)
    clims = [min(data(:)), max(data(:))];
end
if any(isnan(clims)) ||  isequal(clims,[0 0]), clims = [0 1]; end

if ~exist('ttl','var') || isempty(ttl)
    ttl = 'Brainstorm Mesh';
end

% Plot mean amplitude across epochs
if ~exist('fH', 'var') || isempty(fH)
    fH = figure; set(fH, 'Color', 'w', 'Position', [163 483 891 554])
end

% Set up mesh
tH = trimesh(bs_pial_low.Faces,bs_pial_low.Vertices(:,1),bs_pial_low.Vertices(:,2),bs_pial_low.Vertices(:,3));
axis equal; hold on

% Preallocate space for colors to plot (nr vertices x 3 for RGB)
colors = NaN(size(bs_pial_low.Curvature,1),3);
sz     = length(cmap)-1;

% Get curvature
curv = bs_pial_low.Curvature; 

% Implement colors in curvature
colors(curv<=0,:) = .25;
colors(curv>0,:) = .75;

% Get index for data above the requested thresh (default = 0) and select
% those data
ii = find(data>thresh);
Z = data(ii);

% Convert to 1-256
Z_ind = round(sz.*((Z-clims(1)) ./ (clims(2)-clims(1))))+1;

% overlay in colors variable
colors(ii,:) = cmap(Z_ind,:);

% set source pediction as colors
set(tH, 'LineStyle', 'none', 'FaceColor', 'interp', 'FaceVertexCData', colors);
% 
colormap(cmap); colorbar; set(gca, 'CLim',clims, 'view', [-90 0]);

pos = [-.1 0 .1];
light('Position',pos,'Style','local')
% lighting gouraud
material shiny; %dull
title(sprintf('%s', ttl)); 


return