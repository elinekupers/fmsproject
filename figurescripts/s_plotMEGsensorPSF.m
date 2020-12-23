%% s_plotMEGsensorPSF

% subject to plot
subject = 'wlsubj070';

% sensor to plot
sensor = 1;

% Path to brainstorm database and project name
bsDB            = '/Volumes/server/Projects/MEG/brainstorm_db/';
projectName     = 'SSMEG';

% Define vector that can truncate number of sensors 
keep_sensors = logical([ones(157,1); zeros(192-157,1)]); % TODO: Figure out a more generic way to define keep_sensors

% Find Brainstorm anatomy file
d = dir(fullfile(bsDB, projectName, 'data', subject, 'R*'));
bsData = fullfile(d(1).folder, d(1).name);    
bsAnat = fullfile(bsDB, projectName, 'anat', subject);
    
% Load relevant matrices
G_constrained = getGainMatrix(bsData, keep_sensors);

% reduce to sensor of interest
colors = G_constrained(sensor,:)';

% normalize
colors = colors./max(colors);

% shift to be positive
colors = colors+abs(min(colors));

% visualize
clim = [];
thresh = [];
meshType = 'smooth';
ttl = [];
visualizeBrainstormMesh(bsAnat,colors,thresh,clim, meshType,ttl)