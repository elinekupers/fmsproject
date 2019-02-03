%% s_headmodelGeometryDemo.m

subject         = 'wlsubj070';
bsDB            = '/Volumes/server/Projects/MEG/brainstorm_db/';
projectName     = 'SSMEG';
keep_sensors    = logical([ones(157,1); zeros(192-157,1)]); % Note: Figure out a more generic way to define keep_sensors

   
d = dir(fullfile(bsDB, projectName, 'data', subject, 'R*'));
bsData = fullfile(d(1).folder, d(1).name);    
bsAnat = fullfile(bsDB, projectName, 'anat', subject);

fsAnat = fullfile(sprintf('/Volumes/server/Freesurfer_subjects/%s/surf/',subject));

%% 1. Load relevant matrices

G_constrained = getGainMatrix(bsData, keep_sensors);


% Load V1-3 template with unitized phases in downsampled brainstorm format (computed by interp_retinotopy.m)
areas    = load(fullfile(bsAnat, 'areas_overlay.mat')); % [1xNsources] Every value between [-3,3] is inside V1-3, zeros refer to outside of visual cortex. CHECK: Positive values represent LH (?) Negative values RH (?)
eccen    = load(fullfile(bsAnat, 'eccen_overlay.mat')); % [1xNsources] Nonzero value represents vertex preferred eccentricity in degrees, zeros refer to outside of visual cortex
polarang = load(fullfile(bsAnat, 'angle_overlay.mat')); % [1xNsources] Nonzero value represents vertex preferred polar angle in degrees, zeros refer to outside of visual cortex

% Get V1 template
template.V1       = abs(areas.sub_bs_areas)==1;

% For reference, get polar angle of each vertex in degrees (not used yet)
polarang.sub_bs_angle_rad = pi/180*(90-polarang.sub_bs_angle);


paV1 = polarang.sub_bs_angle_rad(template.V1);
ecV1 = eccen.sub_bs_eccen(template.V1);

eccen = 4;
pa    = [1,1.4];

[val,idx] = find((polarang.sub_bs_angle>pa(1))&(polarang.sub_bs_angle<1.4));

thisVertex  = idx(1);



figure; megPlotMap(G_constrained(:,thisVertex),[],[],'bipolar');