% s_plotTemplateBrainstormMesh

% define subject
subject = 'wlsubj039';

% define brainstorm anatomy dir
bsDB        = '/Volumes/server/Projects/MEG/brainstorm_db';
projectName = 'SSMEG';
anatDir     = fullfile(bsDB, projectName, 'anat', subject);

% define brainstorm data dir
% d = dir(fullfile(bsDB, projectName, 'data', subject, 'R*'));
% dataDir = fullfile(d(1).folder, d(1).name);   

area = 'V123';
stimEccen = 11;
template = getTemplate(anatDir, area, 11);

colors = template.(area);
thresh = [];
clims  = [];
meshType = [];
ttl = [];

visualizeBrainstormMesh(anatDir, double(colors)', thresh, clims, meshType, ttl)
