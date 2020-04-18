function [template, polarang, eccen] = getTemplate(anatDir, whichVisualAreas, stimEccen)
% Function to load Brainstorm mesh templates of visual areas. Templates
% should be created from Benson Docker for a particular Freesurfer subject,
% and downsampled to Brainstorm mesh and saved in anat folder in Brainstorm
% database

% INPUTS:
% anatDir            : [str] path to subject's anatomy folder in Brainstorm database
% whichVisualAreas   : [str] string to define whether you want templates
%                            from all early visual areas ['V123'], 'V1',
%                            'V2', or 'V3', all 12 areas in the Benson et 
%                            al. 2018 atlas: 'benson18atlas', or all in
%                            the Wang et al. 2015 atlas: 'wang15atlas'.
%                            Note that the wang atlas does not contain
%                            polar angle/eccen maps, so we will use all
%                            vertices.
% stimEccen          : [int] or [vect] eccentricity in degrees to limit vertices in template.

% OUTPUTS:
% template           :  [struct] contrains downsampled brainstorm meshes with
%                       ones for included vertices and zeros for excluded
%                       vertices
% polarang           :  [struct] contrains downsampled brainstorm meshes with
%                       V1-V3 polar angle map
% eccen              :  [struct] contrains downsampled brainstorm meshes with
%                       V1-V3 eccentricity map

% NB: Brainstorm GUI has to be open
if ~exist('whichVisualAreas','var') || isempty(whichVisualAreas)
    whichVisualAreas = 'V123';
end

if ~exist('stimEccen','var') || isempty(stimEccen)
    stimEccen = [];
end


%%
% [ OLD Template ]
% Load V1-3 template with unitized phases in downsampled brainstorm format (computed by interp_retinotopy.m)
% areas    = load(fullfile(anatDir, 'areas_overlay.mat')); % [1xNsources] Every value between [-3,3] is inside V1-3, zeros refer to outside of visual cortex. Positive values represent dorsal, negative values ventral
% eccen    = load(fullfile(anatDir, 'eccen_overlay.mat')); % [1xNsources] Nonzero value represents vertex preferred eccentricity in degrees, zeros refer to outside of visual cortex
% polarang = load(fullfile(anatDir, 'angle_overlay.mat')); % [1xNsources] Nonzero value represents vertex preferred polar angle in degrees, zeros refer to outside of visual cortex

benson.areas    = load(fullfile(anatDir, 'benson14areas_overlay.mat')); % [1xNsources] Every value between [1 and 3] is inside V1-3, zeros refer to outside of visual cortex, 4-12 are other visual areas
benson.eccen    = load(fullfile(anatDir, 'benson14eccen_overlay.mat')); % [1xNsources] Nonzero value represents vertex preferred eccentricity in degrees, zeros refer to outside of visual cortex, 
benson.polarang = load(fullfile(anatDir, 'benson14angle_overlay.mat')); % [1xNsources] Nonzero value represents vertex preferred polar angle in degrees, zeros refer to outside of visual cortex

wang15.areas      = load(fullfile(anatDir, 'wang15areas_overlay.mat')); 

switch whichVisualAreas
    case 'V123'
        % Get vertices in V1-V3
        template.(whichVisualAreas)     = (benson.areas.sub_bs_areas>0 & benson.areas.sub_bs_areas <=3);
    
    case 'V1'
        % Get only vertices in V1
        template.(whichVisualAreas)     = benson.areas.sub_bs_areas==1; 
        
    case 'V2'
         % Get only vertices in V2
        template.(whichVisualAreas)     = benson.areas.sub_bs_areas==2;
        
    case 'V3'
        % Get only vertices in V3
        template.(whichVisualAreas)     = benson.areas.sub_bs_areas==3;
    
    case 'benson18atlas'
        % Get vertices in all 12 visual areas (1:  'V1', 2: 'V2',  3: 'V3',  4: 'hV4',  5: 'VO1', 6:  'VO2',  7: 'LO1', 8: 'LO2', 9: 'TO1', 10: 'TO2', 11: 'V3b', 12: 'V3a'
        template.(whichVisualAreas)     = benson.areas.sub_bs_areas > 0;
        
    case 'wang15atlas'
        template.(whichVisualAreas)     = wang15.areas.sub_bs_wang > 0;

end

% For reference, get polar angle of each vertex in degrees (not used for now)
benson.polarang.sub_bs_angle_rad = pi/180*(90-benson.polarang.sub_bs_angle);

% define other maps for output
polarang = benson.polarang;
eccen    = benson.eccen;

% Limit to stimulus eccentricity
if ~isempty(stimEccen)
    
    % Add new data to new field with extension    
    tmpName = 'StimEccen';
    
    % We don't have eccentricity values for the wang atlas, so we use all
    % vertices within the ROIs.
    if strcmp(whichVisualAreas, 'wang15atlas')
        eccenMask = ones(size(template.(whichVisualAreas)));
        warning('(%): stimEccen will not be used with wang15atlas', mfilename)
    else
        eccenMask = zeros(size(benson.eccen.sub_bs_eccen));
    
        if length(stimEccen)==1     
            eccenMask((benson.eccen.sub_bs_eccen<=stimEccen(1))) = 1;
        elseif length(stimEccen)==2
            eccenMask((benson.eccen.sub_bs_eccen>=stimEccen(1)) & (benson.eccen.sub_bs_eccen<=stimEccen(2))) = 1;        
        end
        
    end
    
    % Mask the template with eccentricity mask
    template.([whichVisualAreas '_' tmpName]) = template.(whichVisualAreas).*eccenMask;

%     % Debug figure to check number of vertices selected
%     if  strcmp(whichVisualAreas, 'wang15atlas')
%         figure, histogram(wang15.areas.sub_bs_wang(wang15.areas.sub_bs_wang > 0), 'NumBins', 90, 'Normalization', 'cumcount'); hold on;
%         legend({'All vertices, all eccen'});
%     else
%         figure, histogram(eccen.sub_bs_eccen(eccen.sub_bs_eccen>0), 'NumBins', 90, 'Normalization', 'cumcount'); hold on;
%         histogram(eccen.sub_bs_eccen(logical(eccenMask)==1), 'NumBins', 11, 'Normalization', 'cumcount');
%         histogram(eccen.sub_bs_eccen(logical(template.([whichVisualAreas '_' tmpName]))==1), 'NumBins',11);
%         legend({'All vertices, all eccen', 'All vertices, within requested eccen', 'Vertices within visual area, within requested eccen'});
%     end
end

return