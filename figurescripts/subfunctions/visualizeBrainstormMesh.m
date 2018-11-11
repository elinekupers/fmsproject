function visualizeBrainstormMesh(anatDir, colors, thresh, clims, meshType, ttl)



%% Check inputs

if ~exist('meshType','var') || isempty(meshType)
    bs_pial_low = load(fullfile(anatDir, 'tess_cortex_pial_low.mat'));
else
    if strcmp(meshType,'smooth')
        bs_pial_low = load(fullfile(anatDir, 'tess_cortex_pial_low_fig.mat'));
    elseif strcmp(meshType,'unsmooth')
        bs_pial_low = load(fullfile(anatDir, 'tess_cortex_pial_low.mat'));
    end
end

if ~exist('thresh','var') || isempty(thresh)
    thresh = 0;
end

if ~exist('clims','var') || isempty(clims)
    clims = [-1 1].*max(colors(:));
end

if ~exist('ttl','var') || isempty(ttl)
    ttl = 'Brainstorm Mesh';
end

% Define colorbar colors
cmap = [gray(128); jet(128)];

% Plot mean amplitude across epochs
figure; set(gcf, 'Color', 'w', 'Position', [163 483 891 554])

% Set up mesh
tH = trimesh(bs_pial_low.Faces,bs_pial_low.Vertices(:,1),bs_pial_low.Vertices(:,2),bs_pial_low.Vertices(:,3));
axis equal; hold on

curv = bs_pial_low.Curvature - max(bs_pial_low.Curvature);
colors(colors<=thresh) = curv(colors<=thresh);

for tt = 1:size(colors,2)

    % set source pediction as colors
    set(tH, 'LineStyle', 'none', 'FaceColor', 'interp', 'FaceVertexCData',double(colors(:,tt)));

    colormap(cmap); colorbar; set(gca, 'CLim',clims);

    pos = [-.1 0 .1];
    light('Position',pos,'Style','local')
    lighting gouraud
    material shiny; %dull
    title(sprintf('%s, timepoint %d', ttl, tt)); 

    drawnow;
    pause(.1);
end


% Other ways to set curvature colors
% colors = zeros(size(bs_pial_low.Vertices,1),1);
% colors(bs_pial_low.Curvature<0) = -1.5;
% colors(bs_pial_low.Curvature>=0) = -.5;

return