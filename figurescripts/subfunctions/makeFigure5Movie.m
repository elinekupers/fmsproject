function makeFigure5Movie(varargin)
% This is a function to make a movie of simulated synchronous and
% asynchronous cortical time series, propagated to MEG sensors.
% simulation amplitudes are presented in Figure 5 from the manuscript:
%   A visual encoding model links magnetoencephalography
%   signals to neural synchrony in human cortex.
%       by Kupers, Benson, Winawer (YEAR) JOURNAL.
%
%
% INPUTS:
%   [subjectToPlot]        :  (int)  subject nr to plot, default is 12
%   [saveFig]               :  (bool) true/false save figures
%
% Example 1:
%  makeFigure5Movie('subjectToPlot', 1)
%
% By Eline Kupers (NYU) 2017

%% 0. Set up paths and define parameters

p = inputParser;
p.KeepUnmatched = true;
p.addParameter('subjectToPlot', 12);
p.addParameter('headmodelType', 'OS', @(x) any(validatestring(x,{'OS', 'BEM'})));
p.addParameter('highResSurf', false, @islogical);
p.addParameter('area', 'V123', @(x) any(validatestring(x,{'V1', 'V2', 'V3','V23','V123', 'benson18atlas', 'wang15atlas'})));
p.addParameter('eccenLimitDeg', [0.18 11], @isnumeric);
p.addParameter('useConstrainedDipoles', true, @islogical);
p.addParameter('saveMovie',false);
p.parse(varargin{:});


% Rename variables
subjectToPlot         = p.Results.subjectToPlot;
headmodelType         = p.Results.headmodelType;
highResSurf           = p.Results.highResSurf;
area                  = p.Results.area;
eccenLimitDeg         = p.Results.eccenLimitDeg;
useConstrainedDipoles = p.Results.useConstrainedDipoles;
saveMovie             = p.Results.saveMovie;


%% 1. Define subjects and synchrony variables
% Get subject names and corresponding data session number
subject = getSubjectIDs;

% Where to save movie?
if saveMovie
    figureDir = fullfile(fmsRootPath,'figures',subject{subjectsToPlot});
end

% Path to brainstorm database and project name
bsDB            = '/Volumes/server/Projects/MEG/brainstorm_db/';
projectName     = 'SSMEG';

% Number of iterations for the random coherence prediction of the forward model
n        	= 10;        % number of timepoints (ms)
nrEpochs    = 1000;      % number of epochs
theta       = 0;         % von mises mean, equal for three distributions (syn, asyn and mix)
kappa.syn   = 100*pi;    % very narrow von Mises
kappa.asyn  = 0;         % very broad (uniform) von Mises
kappa.mix   = 0.27*pi;   % in-between width size von Mises (note: we are not plotting these values for this figure)

% Define vector that can truncate number of sensors
keep_sensors = logical([ones(157,1); zeros(192-157,1)]); % TODO: Figure out a more generic way to define keep_sensors
    
% Get brainstorm data
d = dir(fullfile(bsDB, projectName, 'data', subject{subjectToPlot}, 'R*'));
bsData = fullfile(d(1).folder, d(1).name);
bsAnat = fullfile(bsDB, projectName, 'anat', subject{subjectToPlot}, 'lowres');

%% 1. Load relevant matrices
G = getGainMatrix(bsData, keep_sensors, headmodelType, highResSurf, useConstrainedDipoles);

% Get V123 template limited to 11 degrees eccentricity
template = getTemplate(bsAnat, area, eccenLimitDeg);

% Simulate coherent, in between or mixture, adn incoherent source time
% series and compute predictions from forward model (w)
tmp = getForwardModelPredictions(G, template.([area '_StimEccen']), [], n, nrEpochs, theta, kappa);

% Define time points to plot
onecyclems = (1/2) * n;
nrCyclesToPlot = 2;
timePoints = 1:1:ceil(nrCyclesToPlot*onecyclems);

% Get mean across epochs for several time points
for t = timePoints
    predictedTS(t,:) = squeeze(mean(tmp.c(:,t,:),3,'omitnan'));
end    


%% Visualize predictions
% Define image size
scr_sz = get(0,'ScreenSize');
rows = ceil(scr_sz(3)/3);
cols = ceil(scr_sz(3)/3);

fH = figure(1); clf; set(gcf, 'Position', [1 1 rows, cols], 'Color', 'w');

% Set colormap, colormap limits
cmap     = bipolar;
clim     = [-1,1].*abs(min(predictedTS(:)));
fprintf('\n(%s): Creating frames.', mfilename);

F(length(timePoints)) = struct('cdata',[], 'colormap', []);
for t = 1:length(timePoints)
    cla;
    [~,ch] = megPlotMap(predictedTS(timePoints(t),:),clim,[],cmap,sprintf('t = %3.0f ms',timePoints(t)),[],[], ...
        'mask', 'convex', ...
        'pointsize', 10);
    set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);
    drawnow;
    
    F(t) = getframe(fH);
end
    
 
% Create movie
mov = implay(F, 4);

if saveMovie
    figureDir = fullfile(fmsRootPath,'figures',subject{subjectToPlot});
    fname = fullfile(figureDir, 'timeseriesForwardModelPredictionMovieFigure5');
    vid = VideoWriter(fname, 'MPEG-4');
    open(vid);
    
    writeVideo(vid, F);
end

fprintf('\n(%s): Done!', mfilename)


