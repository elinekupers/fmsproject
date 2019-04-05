% s_plotTemplateBrainstormMesh

% define subject
subject = 'wlsubj005';

% define brainstorm anatomy dir
bsDB        = '/Volumes/server/Projects/MEG/brainstorm_db';
projectName = 'SSMEG';
anatDir     = fullfile(bsDB, projectName, 'anat', subject);

% What visual area?
area = 'V123';

% What deg of visual angle to limit template?
stimEccen = 11;

% Get template
template = getTemplate(anatDir, area, 11);

% Define plotting params 
colors = template.(area);
thresh = [];
clims  = [];
meshType = [];
ttl = sprintf('%s: Brainstorm mesh, %s', subject, area);

% Plot it!
visualizeBrainstormMesh(anatDir, double(colors)', thresh, clims, meshType, ttl)

% Save the figure
figureDir = fullfile(fmsRootPath, 'figures', subject);
fig_ttl   = sprintf('brainstormMesh_%s_%s', subject, area);
figurewrite(fullfile(figureDir,fig_ttl),[],[1 300],'.',1)