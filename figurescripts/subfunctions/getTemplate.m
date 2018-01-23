function template = getTemplate(anatDir, whichVisualAreas, stimEccen)

% Function to load Brainstorm mesh templates of visual areas. Templates
% should be created from Benson Docker for a particular Freesurfer subject,
% and downsampled to Brainstorm mesh and saved in anat folder in Brainstorm
% database

% INPUTS:
% anatDir            : [str] path to subject's anatomy folder in Brainstorm database
% whichVisualAreas   : [str] string to define whether you want templates
%                            from all visual areas [V1-V3] or just V1
% stimEccen          : [int] eccentricity in degrees to limit vertices in template

% OUTPUTS:
% template           :  [struct] contrains downsampled brainstorm meshes with
%                       ones for included vertices and zeros for excluded
%                       vertices

% NB: Brainstorm GUI has to be open
if ~exist('whichVisualAreas','var') || isempty(whichVisualAreas)
    whichVisualAreas = 'all';
end

if ~exist('stimEccen','var') || isempty(stimEccen)
    stimEccen = [];
end

%%

% Load V1-3 template with unitized phases in downsampled brainstorm format (computed by interp_retinotopy.m)
areas    = load(fullfile(anatDir, 'areas_overlay.mat')); % [1xNsources] Every value between [-3,3] is inside V1-3, zeros refer to outside of visual cortex. CHECK: Positive values represent LH (?) Negative values RH (?)
eccen    = load(fullfile(anatDir, 'eccen_overlay.mat')); % [1xNsources] Nonzero value represents vertex preferred eccentricity in degrees, zeros refer to outside of visual cortex
polarang = load(fullfile(anatDir, 'angle_overlay.mat')); % [1xNsources] Nonzero value represents vertex preferred polar angle in degrees, zeros refer to outside of visual cortex

switch whichVisualAreas
    case 'all'
        % Get only vertices in V1
        template.V123     = abs(areas.sub_bs_areas)>0;
        template.V1       = abs(areas.sub_bs_areas)==1;
        template.V2       = abs(areas.sub_bs_areas)==2;
        template.V3       = abs(areas.sub_bs_areas)==3;
    
    case 'V1'
        % Get only vertices in V1
        template.V1     = abs(areas.sub_bs_areas)==1;     
end

% For reference, get polar angle of each vertex in degrees (not used yet)
polarang.sub_bs_angle_rad = pi/180*(90-polarang.sub_bs_angle);

% Limit to stimulus eccentricity
if ~isempty(stimEccen) 
    
    % Get original names
    names = fieldnames(template); 
   
    % For each name, add new data to new field with extension
    for ii = 1:length(names)        
        tmpName = strcat(names{ii},'StimEccen');      
        template.(tmpName) = template.(names{ii}).*(eccen.sub_bs_eccen<=stimEccen);
    end
    
end

return