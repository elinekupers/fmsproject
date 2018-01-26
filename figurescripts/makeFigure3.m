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
exampleSubject  = 1;

% What's the plotting range for individual example and average across
% subjects?
contourLim      = 0.75; % draw contour line at what fraction of the colormap?
colorbarLim     = 97.5; % percentile of data to use for max/min limits of colorbar

% Number of iterations for the random coherence prediction of the forward
% model
n        = 10;         % number of timepoints (ms)
nrEpochs = 100;        % number of epochs

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
   
    % Take mean amplitude across epochs
    amps.c = abs(fft(tmp.c,[],2));
    amps.i = abs(fft(tmp.i,[],2));
    
    w.V1c(s,:) = mean(amps.c(:,2,:),3);
    w.V1i(s,:) = mean(amps.i(:,2,:),3);
    
end

%% Take mean across subjects

w.V1c_mn = mean(w.V1c,1);
w.V1i_mn = mean(w.V1i,1);

%% Visualize predictions

dataAll      = {w.V1c(exampleSubject,:), w.V1i(exampleSubject,:), w.V1c_mn, w.V1i_mn};
colorMarkers = {'r','b', 'r', 'b'};
ttl          = {'Coherent phase S1', ...
                'Incoherent phase S2', ...
                'Coherent phase Average S1-S6', ...
                'Incoherent phase Average S1-S6'};

fH1 = figure(1); clf; set(fH1,'position',[1,600,1400,800], 'Name', 'Figure 3,  V1 model predictions', 'NumberTitle', 'off');
fH2 = figure(2); clf; subplot(1,2,1); megPlotMap(zeros(1,157)); colormap([1 1 1]);
                      subplot(1,2,2); megPlotMap(zeros(1,157)); colormap([1 1 1]);

for ii = 1:length(dataAll)
    
    dataToPlot = dataAll{ii};
    cLims = [-1 1]*prctile(dataToPlot, colorbarLim);
    
    % Plot predictions
    set(0, 'currentfigure', mod(ii,1)+1);    
    subplot(2,2,ii);
    [~,ch] = megPlotMap(dataToPlot,cLims,fH1,'bipolar',[],[],[], ...
        'isolines', contourLim*max(cLims)*[1 1], ...
        'chanindx', dataToPlot > contourLim*max(cLims), ...
        'pointsymbol', '*', ...
        'pointsize', 10);
    
    c = findobj(gca,'Type','Contour'); c.LineWidth = 4;
    pp = findobj(gca,'Marker','*');
    set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12); title(ttl{ii})
    
    % Plot overlap
    set(0, 'currentfigure', mod(ii,1)+2);     
    subplot(1,2,ceil(ii/2)); hold all;   
    contour(c.XData, c.YData, c.ZData, contourLim*max(cLims)*[1 1], 'LineColor',colorMarkers{ii}, 'LineWidth',4);
    scatter(pp(1).XData,pp(1).YData, 150, colorMarkers{ii},'*'); colorbar off;
    
end

if saveFigures
    set(0, 'currentfigure', fH1);
    figurewrite(fullfile(figureDir,'Figure3_predictionV1_oneVSaverage'),[],0,'.',1);
    set(0, 'currentfigure', fH2);
    figurewrite(fullfile(figureDir,'Figure3_overlap'),[],0,'.',1);
    
    % hgexport(gcf,fullfile(figureDir,'Figure2_SLBB_onesubject_2'))
end


