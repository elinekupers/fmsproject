function G = getGainMatrix(data_dir, keep_sensors, headmodelType, highResFlag, useConstrainedDipoles)

% Function to load gain matrix and if requested constrain it to have only 
% one dipole vector (i.e. perpendicular to the surface).

% This headmodel should be computed by Brainstorm and saved in the data_dir
% that point to the data folder of the subject in the Brainstorm database.

% INPUTS:
% data_dir      : [str]  path to subject's data folder in Brainstorm database
% keep_sensors  : [bool] logical vector that keeps requested sensors and truncates the rest in Gain matrix if requested
% type          : [str]  choose from two types of headmodels: overlapping spheres model ('OS') or ('BEM') boundary element model
% highResFlag   : [bool] true/false use of FS size surface mesh

% OUTPUTS:
% G             :  [matrix] (un)constrained Gain matrix from headmodel. Constrained means that the dipoles are forced to be surface normals (thus perpendicular to the surface)

% NB: Brainstorm GUI has to be open

%%

% Define keep_sensors as empty if not used as input
if nargin < 2
    keep_sensors = []; 
elseif nargin < 3
    headmodelType = 'OS';
elseif nargin < 4
    highResFlag = false;
end

    
% Convert headmodelType to string matching mat file from Brainstorm DB
if strcmp(headmodelType, 'OS'), hm_str = 'os_meg';
elseif strcmp(headmodelType, 'BEM'), hm_str = 'openmeeg';
else, error('(%s): headmodel type does not exist/is not available',mfilename); end

% Get headmodel file name
hm_filename = sprintf('headmodel_surf_%s', hm_str);

% Check high resolution surface flag
if highResFlag
    if strcmp(headmodelType, 'openmeeg')
        error('(%s): high resolution surface is not available for BEM headmodel',mfilename)
    end
    hm_filename = [hm_filename '_02'];
end

% Load headmodel
headmodel = load(fullfile(data_dir, [hm_filename '.mat']));

% Keep all sensors in Gain matrix
if isempty(keep_sensors)
    keep_sensors = ones(size(headmodel.Gain,1),1);
end

% Get Gain matrix and truncate to not-nan sensors
G = headmodel.Gain(keep_sensors,:); % [Nsensors x 3*Nvertices]

% Contrained gain matrix
if useConstrainedDipoles
    G = bst_gain_orient(G, headmodel.GridOrient); % [Nsensors x Nsources], equivalent to size BS pial cortex [1x15002]
else
    %% TO DO
%     how to summarize unconstrained gain matrix
end

return