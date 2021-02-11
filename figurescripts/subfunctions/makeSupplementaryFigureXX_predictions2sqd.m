function makeSupplementaryFigureXX_predictions2sqd()

% NOTE: This function is under construction

% This function uses the forward model to create weights for each timepoint
% of the modeled MEG data. Weights are then put into sqd file. The sqd file
% needs to be imported into the database via Brainstorm GUI, under the same
% R* session, so that the aligned between MRI and MEG is identical and does
% not have to be recomputed.

% Once in Brainstorn, an inverse soluation can be made. We use this inverse
% solution to create a model of time varying sources and project this onto
% the brainstorm mesh.

% So far, the entire routine has only been executed for wl_subj010 and wl_subj002.

% Questions:
% - Are the initial headmodel and the headmodel based on the forward prediction
%   the same?
% - When plotting the inverse model on the Brainstorm mesh, should we define
%   the colors as the absolute values of the source model? (Or more general,
%   how should we deal with the complex numbers at the summary/plotting
%   stage?)


%% 0. Set up paths and define parameters

% Brainstorm Database path
bsDB  = '/Volumes/server/Projects/MEG/brainstorm_db/';

% Freesurfer subject path
fsDir = '/Volumes/server/Freesurfer_subjects/';

% Define project name, subject and data/anatomy folders
project_name    = 'SSMEG';
subject         = 'wl_subj010'; % pick 02, 04, 05, 06, 10, 11

% Get directories
d               = dir(fullfile(bsDB, project_name, 'data', subject, 'R*'));
dataDir         = fullfile(d(1).folder, d(1).name);
anatDir         = fullfile(bsDB, project_name, 'anat', subject);
figureDir       = fullfile(fmsRootPath, 'figures'); % Where to save images?

% Parameters for simulated coherent / incoherent signal
nrEpochs        = 100;     % How many epochs?
n               = 10;       % How many timepoints (ms)
f               = 1;        % Frequency for coherent signal (frequency (Hz))

% Booleans for saving files
saveSQD         = false;    % Save forwardmodel into a sqd file?
saveAsFS        = false;    % Save inverse solution as FS mesh?
saveFigures     = true;     % Save figures in the figure folder?

% Define vector that can truncate number of sensors 
keep_sensors = logical([ones(157,1); zeros(192-157,1)]); % Note: Figure out a more generic way to define keep_sensors


%% -----------------------------------------------------------------------
%  1. Get gain matrix and templates for forward model, compute predictions
%  -----------------------------------------------------------------------

G_constrained = getGainMatrix(dataDir, keep_sensors);

% Get V1 template limited to 11 degrees eccentricity
template = getTemplate(anatDir, 'V1', 11);

% Simulate coherent and incoherent source time series and compute
% predictions from forward model (w)
w.V1 = getForwardModelPredictions(G_constrained, template.V1StimEccen, f, n, nrEpochs);

% Summarize by convering sensor time series to amplitudes
tmp = abs(fft(w.V1.c, [], 2));
w.V1.cAmp = tmp(:, f+1, :);

tmp = abs(fft(w.V1.i, [], 2));
w.V1.iAmp = tmp(:, f+1, :);

% Take mean across epochs
w.V1.iAmp_mn = mean(w.V1.iAmp,3);
w.V1.cAmp_mn = mean(w.V1.cAmp,3);


%% -----------------------------------------------------------------------
%  2. Visualize mean forward model predictions
%  -----------------------------------------------------------------------

figure(1),
clims = max([w.V1.cAmp_mn; w.V1.iAmp_mn]) * [-1 1];
subplot(1,2,1); megPlotMap(w.V1.cAmp_mn, clims, [], 'bipolar', 'Coherent phase', [], [], 'isolines', 0.3*max(clims) * [1 1]);
subplot(1,2,2); megPlotMap(w.V1.iAmp_mn, clims, [], 'bipolar', 'Incoherent phase', [], [], 'isolines', 0.15*max(clims) * [1 1]);

 if saveFigures
    hgexport(gcf, fullfile(figureDir, sprintf('SF2A_predictionForwardModelV1wtimeseries_%s.eps',subject)))
 end
 
%% -----------------------------------------------------------------------
%  3. Create SQD file (Only necessary once, after that it should be saved in Brainstorm DB)
%  -----------------------------------------------------------------------

% We create a new SQD file with the prediction from the forward model.
% This new sqd file needs to be imported into Brainstorm session under the
% same subject in data tab to complete the inverse model step (4)

if saveSQD
    
    ok = savePrediction2SQD(fmsRootPath, subject, w.V1.c, 'V1ForwardCoherent');
    ok = savePrediction2SQD(fmsRootPath, subject, w.V1.i, 'V1ForwardIncoherent');

end

%% -----------------------------------------------------------------------
%  4. Create inverse prediction (s) from forward model (w)
%  -----------------------------------------------------------------------

d_inverse = dir(fullfile(dataDir, 'results_MN_MEG_KERNEL*18*'));

% Get inverse model for coherent signals
inverseC  = fullfile(d_inverse(1).folder, d_inverse(1).name);
s.V1.c = getInverse(inverseC, w.V1.c, keep_sensors);

% Get inverse model for coherent signals
inverseI  = fullfile(d_inverse(2).folder, d_inverse(2).name);
s.V1.i = getInverse(inverseI, w.V1.i, keep_sensors);

% Summarize predicted source amplitudes

tmp = abs(fft(s.V1.c, [], 2));
s.V1.cAmp = tmp(:, f+1, :);

tmp = abs(fft(s.V1.i, [], 2));
s.V1.iAmp = tmp(:, f+1, :);

%% -----------------------------------------------------------------------
%  5. Visualize source predictions (from inverse model)
%  -----------------------------------------------------------------------

s_all = struct2cell(s.V1);

% Plotting params
labels   = {'Coherent', 'Incoherent'};
meshType = 'smooth';   
thresh   = 0;

for source = [1,2] % 1 is coherent (all vertices have the same phase of a 12 Hz sine), 2 is incoherent (all vertices have a random phase)
    
    ttl = sprintf('Inverse solution for %s signal, threshold: %1.2f', labels{source}, thresh);
    
    thisSource = s_all{source};
    
    thisSource = reshape(thisSource,size(thisSource,1),[]);
    
    colors = abs(thisSource);
    
    % Show sources for timepoints of 2 epochs
    visualizeBrainstormMesh(anatDir, colors(:,1:(2*n)), thresh, [], meshType, ttl)
    
    % Show mean timeseries   
%     %Plot mean timeseries of V1 sources
%     figure; subplot(2,1,1);
%     plot(t, mean(s_all{2+source}(logical(template.V1StimEccenAmplitudes'),:,:),3)); hold on;
%     plot(t, mean(mean(s_all{2+source}(logical(template.V1StimEccenAmplitudes'),:,:),3)), 'k-', 'LineWidth', 4)
%     xlabel('Time (s)'); ylabel('Source amplitude (??)'); box off; set(gca,'TickDir', 'out')
%     title(sprintf('%s: Mean timeseries of V1 sources', labels{source}))
%     
%     % Plot mean timeseries of non visual sources
%     nonVisualVertices = abs(areas.sub_bs_areas)==0;
%     subplot(2,1,2);
%     plot(t,mean(s_all{2+source}(nonVisualVertices,:,:),3)); hold on;
%     plot(t, mean(mean(s_all{2+source}(nonVisualVertices,:,:),3)), 'k-', 'LineWidth', 4)
%     xlabel('Time (s)'); ylabel('Source amplitude (??)'); box off; set(gca,'TickDir', 'out')
%     title(sprintf('%s: Mean timeseries of non V1 sources', labels{source}))
%     
%     if saveFigures
%         hgexport(gcf, fullfile(figureDir, sprintf('Figure5B_sourcetimeseriesV1_%s_%s.eps', labels{source}, subject)))
%     end
    
end

%% Show mean amplitudes on brainstorm mesh

thresh   = 0.03;
clims    = [-1 1]* max([mean(s.V1.cAmp,3); mean(s.V1.iAmp,3)]);

% --- Coherent ---
ttl      = sprintf('Mean 1Hz amplitude: Coherent, threshold: %1.2f', thresh);
colors   = mean(s.V1.cAmp,3);

visualizeBrainstormMesh(anatDir, colors, thresh, clims, meshType, ttl)

% --- Incoherent ---
ttl      = sprintf('Mean 1Hz amplitude: Incoherent, threshold: %1.2f', thresh);
colors   = mean(s.V1.iAmp,3);

visualizeBrainstormMesh(anatDir, colors, thresh, clims, meshType, ttl)


%% -----------------------------------------------------------------------
%  5. If requested: Save BS mesh to FS space
%  -----------------------------------------------------------------------

if saveAsFS
    
    % NB: BRAINSTORM GUI HAS TO BE OPEN FOR THIS STEP
    
    if ~exist('MRIwrite')
        addpath(genpath('/Applications/freesurfer/matlab'));
        addpath(genpath('/Applications/freesurfer/fsfast/toolbox'));
    end
    
    % Set Brainstorm mesh back to FS mesh
    
    % Coherent
    [fs_lh_overlay, fs_rh_overlay] = tess_bst2fs(subject, fullfile(fsDir, subject), mean(s.V1.cAmp,3));
    
    MRIwrite(struct('vol', fs_lh_overlay), fullfile(fsDir, subject,'surf','lh.inverse_V1_coherent.mgz'));
    MRIwrite(struct('vol', fs_rh_overlay), fullfile(fsDir, subject,'surf','rh.inverse_V1_coherent.mgz'));
    
    % Incoherent
    [fs_lh_overlay, fs_rh_overlay] = tess_bst2fs(subject, fullfile(fsDir, subject), mean(s.V1.iAmp,3));
    
    MRIwrite(struct('vol', fs_lh_overlay), fullfile(fsDir, subject,'surf','lh.inverse_V1_incoherent.mgz'));
    MRIwrite(struct('vol', fs_rh_overlay), fullfile(fsDir, subject,'surf','rh.inverse_V1_incoherent.mgz'));

end

return

