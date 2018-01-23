function makeFigure3()

% This is a function to make Figure 3 from the manuscript about forward
% modeling coherent and incoherent neural sources to MEG responses.

% This figure shows the MEG forward model based on coherent and incoherent
% predictions coming from vertices located in V1.

% To runs this script, you need:     
% (1) Access to the SSMEG folder in the brainstorm data base
% (2) MEG_utils and Fieldtrip toolbox added to the paths. For example:
%     tbUse('ForwardModelSynchrony');
%        or to only add the MEG_utils toolbox:
%     addpath(genpath('~/matlab/git/toolboxes/meg_utils'))


%% 0. Set up paths and define parameters

% Path to brainstorm database
bsDB            = '/Volumes/server/Projects/MEG/brainstorm_db/';
figureDir       = fullfile(fmsRootPath, 'figures'); % Where to save images?
saveFigures     = true;     % Save figures in the figure folder?

% Define project name, subject and data/anatomy folders
projectName    = 'SSMEG';

% Which subjects to average?
subject         = {'wl_subj002','wl_subj004','wl_subj005','wl_subj006','wl_subj010','wl_subj011'};

% What's the plotting range for individual example and average across
% subjects?
climsOne = [-.5 .5]*1E-4;
climsAve = [-.5 .5]*1E-4;

% Number of iterations for the random coherence prediction of the forward
% model
n        = 1000;     % number of timepoints (ms)
nrEpochs = 1;        % number of epochs

% Define vector that can truncate number of sensors 
keep_sensors = logical([ones(157,1); zeros(192-157,1)]); % Note: Figure out a more generic way to define keep_sensors


for s = 1:length(subject)
    
    d = dir(fullfile(bsDB, projectName, 'data', subject{s}, 'R*'));
    dataDir = fullfile(d(1).folder, d(1).name);    
    anatDir = fullfile(bsDB, projectName, 'anat', subject{s});
    
    %% 1. Load relevant matrices
    
    G_constrained = getGainMatrix(dataDir, keep_sensors);

    % Get V1 template limited to 11 degrees eccentricity
    template = getTemplate(anatDir, 'V1', 11);

    % Simulate coherent and incoherent source time series and compute
    % predictions from forward model (w)
    tmp = getForwardModelPredictions(G_constrained, template.V1StimEccen, [], n, nrEpochs);
   
    % Take mean across epochs
    w.V1c(s,:) = mean(abs(tmp.c),2);
    w.V1i(s,:) = mean(abs(tmp.i),2);
    
end

%% 3. Visualize predictions from forward model

% plot for subj 002
fH1 = figure(1); clf; set(fH1, 'Position', [1 1 1600 513], 'Name', 'Figure 3: V1 model predictions');

ch = subplot(221); 
megPlotMap(abs(w.V1c(1,:)),climsOne,[],bipolar,[],[],[],'isolines', 0.75*max(climsOne) * [1 1]);
set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12); title('Coherent phase S1')

ch = subplot(222); 
megPlotMap(abs(w.V1i(1,:)),climsOne,[],bipolar,[],[],[],'isolines',  0.75*max(climsOne) * [1 1])
set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12); title('Incoherent phase S1')

%% Compute the mean across subjects

w.V1c_mn = mean(w.V1c,1);
w.V1i_mn = mean(w.V1i,1);

ch = subplot(223); 
megPlotMap(abs(w.V1c_mn),climsAve,[],bipolar,[],[],[],'isolines', 0.75*max(climsAve) * [1 1])
set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12); title('Coherent phase Average (S1-S6)')

ch = subplot(224); 
megPlotMap(abs(w.V1i_mn),climsAve,[],bipolar,[],[],[],'isolines',  0.75*max(climsAve) * [1 1])
set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12); title('Incoherent phase Average (S1-S6)')

if saveFigures
    hgexport(gcf,fullfile(figureDir, 'Figure3_predictionV1_oneVSaverage.eps'));
%     figurewrite(fullfile(figureDir, ['predictionV1_oneVSaverage']), [],0,'.',1);
end
