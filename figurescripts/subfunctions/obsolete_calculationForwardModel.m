
% Number of iterations for the random coherence prediction of the forward
% model
iter = 1000;

for s = 1:length(subject)
    
    d = dir(fullfile(bs_db, project_name, 'data', subject{s}, 'R*'));
    data_dir = fullfile(d(1).folder, d(1).name);
    anat_dir = fullfile(bs_db, project_name, 'anat', subject{s});
    
    %% 1. Load relevant matrices
    % Load Gain matrix created by brainstorm
    headmodel = load(fullfile(data_dir, 'headmodel_surf_os_meg.mat'));
    G = headmodel.Gain(1:157,:); % Gain matrix, [Nsensors x 3*Nvertices]
    
    % Take absolute values of G_contrained - no cancellation possible
    G_constrained = abs(bst_gain_orient(G, headmodel.GridOrient)); % Contrained gain matrix [Nsensors x Nsources], equivalent to size pial cortex [1x15002]
    
    % Load V1-3 template in downsampled brainstorm format (computed by interp_retinotopy.m)
    areas = load(fullfile(anat_dir, 'areas_overlay.mat')); % [1xNsources] Every value between [-3,3] is inside V1-3, zeros refer to outside of visual cortex
    template.V1 = abs(areas.sub_bs_areas)==1;
    
    %     V123PhaseScrambledTemplate = load(fullfile(anat_dir, sprintf('areas_overlay_%s.mat',iterations))); % [1xNsources] Every value between [-3,3] is inside V1-3, zeros refer to outside of visual cortex
    
    % Load V1-3 eccentricity in downsampled brainstorm format
    eccen    = load(fullfile(anat_dir, 'eccen_overlay.mat')); % [1xNsources] Every value from [1 3] is inside V1-3, zeros refer to outside of visual cortex
    polarang = load(fullfile(anat_dir, 'angle_overlay.mat'));
    polarang.sub_bs_angle_rad = pi/180*(90-polarang.sub_bs_angle);
    
    % Restrict V123 to 11 degrees eccentricity (stimulus diameter = 22)
    template.V1StimEccenAmplitudes = template.V1 .*(eccen.sub_bs_eccen<=11);
    
    template.V1StimEccenPhaseCoherent     = template.V1StimEccenAmplitudes * 0;
    template.V1StimEccenPhaseIncoherent   = template.V1StimEccenAmplitudes .* (rand(iter,size(template.V1StimEccenAmplitudes,2)) * 2*pi);
    
    %% 2. Compute forward solution (Method 1)
    % Compute the sensor weights, w, from V1-V3 using the contrained gain field (forward model)
    %     w = G_constrained*V1template'; %  Nsensors x 1;
    
    phAmp2complex = @(r,th) r .* exp(1i*th);
    
    % Make a complex number with amplitudes and phases:
    template.V1coherent = phAmp2complex(template.V1StimEccenAmplitudes,template.V1StimEccenPhaseCoherent);
    
    template.V1incoherent = phAmp2complex(repmat(template.V1StimEccenAmplitudes,[iter,1]),template.V1StimEccenPhaseIncoherent);
    
    % Compute prediction of forward model for template restricted by stimulus eccentricity
    w.V1c(s,:) = G_constrained*template.V1coherent'; %  Nsensors x 1;
    w.V1i(s,:) = mean(abs(G_constrained*template.V1incoherent'),2); %  Nsensors x 1;
    
end