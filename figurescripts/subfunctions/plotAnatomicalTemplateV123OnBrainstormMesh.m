function plotAnatomicalTemplateV123OnBrainstormMesh(subject, area, stimEccen, saveFig)
%
% Function to plot Benson 14 anatomical template of visual cortex (V1-V3)
% on brainstorm mesh (thus downsampled and hemispheres combined)

% INPUTS
% subject   :   (string) name of subject (for example 'wlsubj005')
% area      :   (string) visual area (for example 'V123')
% stimEccen :   (double or vector) eccentricity limits of stimulus in
%                        visual field in deg. Vector indicates upper and lower
%                        boundary, a double indicates lower boundary = 0 and 
%                        upper boundary is X. (for example [2 4] or [11] to get 0-11)
% saveFig   :   (bool)   save figure in rootpath figure dir
%
% Example 1:
% plotAnatomicalTemplateV123OnBrainstormMesh('wlsubj070', 'V1', [0 1], 1)
%

% define brainstorm anatomy dir
bsDB        = '/Volumes/server/Projects/MEG/brainstorm_db';
projectName = 'SSMEG';
anatDir     = fullfile(bsDB, projectName, 'anat', subject);

if ~exist('saveFig', 'var'); saveFig=0; end

if length(stimEccen)==1
    stimEccen(2) = stimEccen(1);
    stimEccen(1) = 0;
end

% Get template
template = getTemplate(anatDir, area, stimEccen);

% Define plotting params 
colors = double(template.V123)./max(template.V123); %template.([area '_StimEccen']);
thresh = [];
clims  = [];
meshType = 'smooth'; %'unsmooth'
ttl = sprintf('%s: Brainstorm mesh, %s', subject, area);

% Plot it!
visualizeBrainstormMesh(anatDir, colors', thresh, clims, meshType, ttl)

% Save the figure
if saveFig
    figureDir = fullfile(fmsRootPath, 'figures', subject);
    fig_ttl   = sprintf('brainstormMesh_%s_%s_%d-%d', subject, area, stimEccen(1), stimEccen(2));
    figurewrite(fullfile(figureDir,fig_ttl),[],[1 300],'.',1)
end