function G_constrained = getGainMatrix(data_dir, keep_sensors)

% Function to load gain matrix and constrain it to have only one dipole 
% vector (that is perpendicular to the surface).

% This headmodel should be computed by Brainstorm and saved in the data_dir
% that point to the data folder of the subject in the Brainstorm database.

% INPUTS:
% data_dir      : [str] path to subject's data folder in Brainstorm database
% max_size      : [bool] logical that keeps requested sensors and truncates the rest in Gain matrix if requested

% OUTPUTS:
% G_constrained :  contrained Gain matrix from headmodel 

% NB: Brainstorm GUI has to be open

%%

% Define keep_sensors as empty if not used as input
if nargin < 2
    keep_sensors = []; 
end
    
% Load headmodel from Brainstorm
headmodel = load(fullfile(data_dir, 'headmodel_surf_os_meg.mat'));

% Keep all sensors in Gain matrix
if isempty(keep_sensors)
    keep_sensors = ones(size(headmodel.Gain,1),1);
end

% Get Gain matrix and truncate to not-nan sensors
G = headmodel.Gain(keep_sensors,:); % [Nsensors x 3*Nvertices]

% Contrained gain matrix
G_constrained = bst_gain_orient(G, headmodel.GridOrient); % [Nsensors x Nsources], equivalent to size BS pial cortex [1x15002]



return