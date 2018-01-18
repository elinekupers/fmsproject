function makeFigure5()

% This function is under construction

% This function uses the forward model to create weights for each timepoint
% of the modeled MEG data. Weights are then put into sqd file, which needs
% to be imported into the database via Brainstorm GUI. So far, this has
% only been done for wl_subj010.

% Once in Brainstorn, an inverse soluation can be made. We use this inverse
% solution to create a model of time varying sources and project this onto
% the brainstorm mesh.



% Questions:
% - Are the initial headmodel and the headmodel based on the forward prediction
%   the same?
% - When plotting the inverse model on the Brainstorm mesh, should we define
%   the colors as the absolute values of the source model? (Or more general,
%   how should we deal with the complex numbers at the summary/plotting
%   stage?)


%% 0. Set up paths and define parameters

% Brainstorm Database path
bs_db = '/Volumes/server/Projects/MEG/brainstorm_db/';

% Define project name, subject and data/anatomy folders
project_name = 'SSMEG';
subject = 'wl_subj010'; % pick 02, 04, 05, 06, 10, 11
% iterations = 'phase_0'; % iterations for the phase scrambled predictioin (number stands for smoothing iterations - zero = no smoothing)

d = dir(fullfile(bs_db, project_name, 'data', subject, 'R*'));

data_dir = fullfile(d(1).folder, d(1).name);
anat_dir = fullfile(bs_db, project_name, 'anat', subject);

% How many epochs?
nrEpochs = 1000;

% Save forwardmodel into a sqd file?
saveSQD = false;

figureDir       = fullfile(fmsRootPath, 'figures'); % Where to save images?
saveFigures     = true;     % Save figures in the figure folder?

%% 1. Make predictions from forwardmodel

% Load headmodel from Brainstorm
headmodel = load(fullfile(data_dir, 'headmodel_surf_os_meg.mat'));

% Get Gain matrix and truncate to first 157 sensors
G = headmodel.Gain(1:157,:); % [Nsensors x 3*Nvertices]

% Contrained gain matrix
G_constrained = bst_gain_orient(G, headmodel.GridOrient); % [Nsensors x Nsources], equivalent to size BS pial cortex [1x15002]

% Load V1-3 template with unitized phases in downsampled brainstorm format (computed by interp_retinotopy.m)
areas    = load(fullfile(anat_dir, 'areas_overlay.mat')); % [1xNsources] Every value between [-3,3] is inside V1-3, zeros refer to outside of visual cortex. CHECK: Positive values represent LH (?) Negative values RH (?)
eccen    = load(fullfile(anat_dir, 'eccen_overlay.mat')); % [1xNsources] Nonzero value represents vertex preferred eccentricity in degrees, zeros refer to outside of visual cortex
polarang = load(fullfile(anat_dir, 'angle_overlay.mat')); % [1xNsources] Nonzero value represents vertex preferred polar angle in degrees, zeros refer to outside of visual cortex

% Get only vertices in V1
template.V1     = abs(areas.sub_bs_areas)==1;

% Get polar angle of each vertex in degrees
polarang.sub_bs_angle_rad = pi/180*(90-polarang.sub_bs_angle);

% Limit to 11 degrees eccentricity
template.V1StimEccenAmplitudes = template.V1.*(eccen.sub_bs_eccen<=11);

%% Make source time series

% Create two cycles of a sine wave
n  = 24;        % number of time points
t  = (1:n)/n;   % s
f  = 2;         % frequency (Hz)

signalCoherent = zeros([size(template.V1StimEccenAmplitudes,2), n, nrEpochs]);
signalIncoherent = zeros(size(signalCoherent));

for epoch = 1:nrEpochs
    
    % Coherent signal gets a fixed phase for each vertex, in one epoch (but
    % different across epochs)
    thisEpochCoherentPhase = rand*2*pi;
    signalCoherent(:,:,epoch) = template.V1StimEccenAmplitudes' * sin(2*pi*f*t + thisEpochCoherentPhase);
    
    % Incoherent signal gets a random phase every vertex and every epoch
    for ii = find(template.V1StimEccenAmplitudes)
        signalIncoherent(ii,:,epoch) = sin(2*pi*f*t + rand*2*pi);
    end
    
    %% Predicted sensor time series
    w.V1c(:,:,epoch) = G_constrained*signalCoherent(:,:,epoch);
    w.V1i(:,:,epoch) = G_constrained*signalIncoherent(:,:,epoch);
    
end



%% Summarize by convering sensor time series to amplitudes
tmp = abs(fft(w.V1c, [], 2));
w.V1cAmp = tmp(:, f+1, :);

tmp = abs(fft(w.V1i, [], 2));
w.V1iAmp = tmp(:, f+1, :);

% Take mean across epochs
w.V1iAmp_mn = mean(w.V1iAmp,3);
w.V1cAmp_mn = mean(w.V1cAmp,3);

%% Visualize
figure(1),
clims = max([w.V1cAmp_mn; w.V1iAmp_mn]) * [-1 1];
subplot(1,2,1); megPlotMap(w.V1cAmp_mn, clims, [], 'bipolar', 'Coherent phase', [], [], 'isolines', 0.4*max(clims) * [1 1]);
subplot(1,2,2); megPlotMap(w.V1iAmp_mn, clims, [], 'bipolar', 'Incoherent phase', [], [], 'isolines', 0.2*max(clims) * [1 1]);

if saveFigures
    hgexport(gcf, fullfile(figureDir, 'Figure5A_predictionForwardModelV1wtimeseries.eps'))
end

%% 3. Create SQD file (Only necessary once, after that it should be saved in Brainstorm DB)

% We create a new SQD file with the prediction from the forward model.
% This new sqd file needs to be imported into Brainstorm session under the
% same subject in data tab to complete the inverse model step (4)

if saveSQD
    
    % Get example meg sqd file (TODO: Make an example sqd file that is small
    % and easy to download)
    [~, meg_files] = meg_load_sqd_data('/Volumes/server/Projects/MEG/SSMEG/09_SSMEG_06_27_2014_wl_subj010/raw/','V1ForwardStimEccen_12HzCycle');
    
    if ~exist(fullfile(fmsRootPath, 'data', subject),'dir')
        mkdir(fullfile(fmsRootPath, 'data', subject));
    end
    
    % Add additional NaNs to pad the sensor space to 192
    dataToSave = reshape(w.V1c, size(w.V1c,1), []);
    dataToSave = [dataToSave; NaN(192-size(dataToSave,1),size(dataToSave,2))];
    
    % For coherent phase (12 Hz)
    newFile = fullfile(fmsRootPath, 'data', subject, 'V1ForwardCoherent.sqd');
    sqdwrite(fullfile(meg_files.folder,meg_files.name),newFile, dataToSave');
    
    % Add additional NaNs to pad the sensor space to 192
    dataToSave = reshape(w.V1i, size(w.V1i,1), []);
    dataToSave = [dataToSave; NaN(192-size(dataToSave,1),size(dataToSave,2))];
    
    % For incoherent phase (random)
    newFile = fullfile(fmsRootPath, 'data', subject, 'V1ForwardIncoherent.sqd');
    sqdwrite(fullfile(meg_files.folder,meg_files.name),newFile,dataToSave');
    
end

%% 4. Create inverse model

% Load Brainstorm downsampled pial surface
bs_pial_low = load(fullfile(anat_dir, 'tess_cortex_pial_low.mat'));

d_inverse_c = dir(fullfile(bs_db, project_name, 'data', subject, '*Coherent*', 'results_MN_MEG_KERNEL*'));
d_inverse_i = dir(fullfile(bs_db, project_name, 'data', subject, '*Incoherent*', 'results_MN_MEG_KERNEL*'));

% Load Brainstorm inverse model
inverseC  = load(fullfile(d_inverse_c.folder, d_inverse_c.name));
inverseI  = load(fullfile(d_inverse_i.folder, d_inverse_i.name));

w.V1c_timeseries = reshape(w.V1c, size(w.V1c,1), []);
w.V1i_timeseries = reshape(w.V1i, size(w.V1i,1), []);


% Create source response for each timepoint
for ii = 1:size(w.V1c_timeseries,2)
    s.V1c_timeseries(:,ii) = inverseC.ImagingKernel*w.V1c_timeseries(:,ii);
    s.V1i_timeseries(:,ii) = inverseI.ImagingKernel*w.V1i(1:157,ii);
end

% Reshape back into epoched source time series
s.V1c = reshape(s.V1c_timeseries, [size(s.V1c_timeseries,1), size(w.V1c,2), size(w.V1c,3)]);
s.V1i = reshape(s.V1i_timeseries, [size(s.V1i_timeseries,1), size(w.V1i,2), size(w.V1i,3)]);


%% Summarize predicted source amplitudes

tmp = abs(fft(s.V1c, [], 2));
s.V1cAmp = tmp(:, f+1, :);

tmp = abs(fft(s.V1i, [], 2));
s.V1iAmp = tmp(:, f+1, :);

%% Visualize
s_all = struct2cell(s);

labels = {'Coherent', 'Incoherent'};

for source = [1,2] % 1 is coherent (all vertices have the same phase of a 12 Hz sine), 2 is incoherent (all vertices have a random phase)
    
    thisSource = s_all{source};
    
    %% Show sources for each time point
    
    % Visualize: set up curvature colors
    colors = zeros(size(bs_pial_low.Vertices,1),1);
    colors(bs_pial_low.Curvature<0) = -1.5;
    colors(bs_pial_low.Curvature>=0) = -.5;
    
    % Visualize: define colorbar colors
    cmap = [gray(128); jet(128)];
    
    % Visualize: set up mesh
    figure; tH = trimesh(bs_pial_low.Faces,bs_pial_low.Vertices(:,1),bs_pial_low.Vertices(:,2),bs_pial_low.Vertices(:,3));
    axis equal; hold on
    
    % Plot a subset of timepoints
    for ii = 1:(size(thisSource,2)/nrEpochs)*5
        
        % Define colors as the absolute values of the source model
        colors = abs(thisSource(:,ii));
        
        % set source pediction as colors
        set(tH, 'LineStyle', 'none', 'FaceColor', 'interp', 'FaceVertexCData',colors);
        drawnow;
        colormap(cmap); set(gca, 'CLim',2*pi*[-.01 .01]); colorbar;
        
        pos = [-.1 0 .1];
        %         light('Position',pos,'Style','local')
        lighting gouraud
        material shiny; %dull
        title(ii);
        pause(.1);
    end
    
    
    %% Show mean timeseries
    
    %Plot mean timeseries of V1 sources
    figure; subplot(2,1,1);
    plot(t, mean(s_all{2+source}(logical(template.V1StimEccenAmplitudes'),:,:),3)); hold on;
    plot(t, mean(mean(s_all{2+source}(logical(template.V1StimEccenAmplitudes'),:,:),3)), 'k-', 'LineWidth', 4)
    xlabel('Time (s)'); ylabel('Source amplitude (??)'); box off; set(gca,'TickDir', 'out')
    title(sprintf('%s: Mean timeseries of V1 sources', labels{source}))
    
    % Plot mean timeseries of non visual sources
    nonVisualVertices = abs(areas.sub_bs_areas)==0;
    subplot(2,1,2);
    plot(t,mean(s_all{2+source}(nonVisualVertices,:,:),3)); hold on;
    plot(t, mean(mean(s_all{2+source}(nonVisualVertices,:,:),3)), 'k-', 'LineWidth', 4)
    xlabel('Time (s)'); ylabel('Source amplitude (??)'); box off; set(gca,'TickDir', 'out')
    title(sprintf('%s: Mean timeseries of non V1 sources', labels{source}))
    
    if saveFigures
        hgexport(gcf, fullfile(figureDir, sprintf('Figure5B_sourcetimeseriesV1_%s.eps', labels{source})))
    end
    
end

%% Show mean amplitudes

thresh = 0.01;

% Plot mean amplitude across epochs
figure; set(gcf, 'Color', 'w', 'Position', [1 484 2407 554])

% ----- Coherent ----
subplot(1,2,1);
tH = trimesh(bs_pial_low.Faces,bs_pial_low.Vertices(:,1),bs_pial_low.Vertices(:,2),bs_pial_low.Vertices(:,3));
axis equal; hold on

% Define colors as the absolute values of the source model
colors = abs(mean(s.V1cAmp,3));
colors(colors<=thresh) = -.01;

% set source pediction as colors
set(tH, 'LineStyle', 'none', 'FaceColor', 'interp', 'FaceVertexCData',colors);
drawnow;
colormap(cmap); set(gca, 'CLim',2*pi*[-.01 .01]); colorbar;

pos = [-.1 0 .1];
light('Position',pos,'Style','local')
lighting gouraud
material shiny; %dull
title('Mean 2Hz amplitude: Coherent');

% ----- Incoherent ----
subplot(1,2,2);
tH = trimesh(bs_pial_low.Faces,bs_pial_low.Vertices(:,1),bs_pial_low.Vertices(:,2),bs_pial_low.Vertices(:,3));
axis equal; hold on

% Define colors as the absolute values of the source model
colors = abs(mean(s.V1iAmp,3));
colors(colors<=thresh) = -.01;

% set source pediction as colors
set(tH, 'LineStyle', 'none', 'FaceColor', 'interp', 'FaceVertexCData',colors);
drawnow;
colormap(cmap); set(gca, 'CLim',2*pi*[-.01 .01]); colorbar;

pos = [-.1 0 .1];
light('Position',pos,'Style','local')
lighting gouraud
material shiny; %dull
title('Mean 2Hz amplitude: Incoherent');

if saveFigures
    hgexport(gcf, fullfile(figureDir, 'Figure5C_InverseSolutionForV1.eps'))
end

return
%% 5. NOT READY YET: Save to FS space

% NB: BRAINSTORM GUI HAS TO BE OPEN FOR THIS STEP

% if ~exist('MRIwrite')
%     addpath(genpath('/Applications/freesurfer/matlab'));
%     addpath(genpath('/Applications/freesurfer/fsfast/toolbox'));
% end
%
%
% % Where are the FS subjects?
% fs_dir = '/Volumes/server/Freesurfer_subjects/';
%
% % Set Brainstorm mesh back to FS mesh
% [fs_lh_overlay fs_rh_overlay] = tess_bst2fs(subject, fullfile(fs_dir, subject), s.V1c);
%
%
% MRIwrite(struct('vol', fs_lh_overlay), fullfile(fs_dir, subject,'surf','lh.inverse_V1_coherent.mgz'));
% MRIwrite(struct('vol', fs_rh_overlay), fullfile(fs_dir, subject,'surf','rh.inverse_V1_coherent.mgz'));

