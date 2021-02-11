function makeFigure3Movie(varargin)
% This is a function to make a movie of the stimulus-locked and broadband
% data presented in Figure 3 from the manuscript:
%   A visual encoding model links magnetoencephalography
%   signals to neural synchrony in human cortex.
%       by Kupers, Benson, Winawer (YEAR) JOURNAL.
%
%
% INPUTS:
%   [subjectsToPlot]        :  (int)  subject nr to plot, default is 12
%   [saveFig]               :  (bool) true/false save figures
%
% Example 1:
%  makeFigure3Movie('subjectsToPlot', 1, 'saveFig', true)
%
% By Eline Kupers (NYU) 2017

%% 0. Set up paths and define parameters

p = inputParser;
p.KeepUnmatched = true;
p.addParameter('subjectsToPlot', 2);
p.addParameter('saveFig', true, @islogical);
p.parse(varargin{:});

% Rename variables
subjectsToPlot        = p.Results.subjectsToPlot;
saveFig               = p.Results.saveFig;

% Where do data live?
dataDir = fullfile(fmsRootPath, 'data'); 

%% 1. Define subjects and synchrony variables

% Get subject names and corresponding data session number
[subject, dataSession] = getSubjectIDs;
whichSession = dataSession(subjectsToPlot);

% Load data
data     = loadData(fullfile(dataDir, subject{subjectsToPlot}), whichSession, 'type', 'timeseries');
dataFull = data.ts(:,:,data.condEpochsFull==1);
dataBlank = data.ts(:,:,data.condEpochsBlank==1);

% Average time series across full field epochs
tsFullAverage = mean(dataFull,3);
tsBlankAverage = mean(dataBlank,3);

% Add bad channels back in as nans
tsFullAverage = to157chan(tsFullAverage', ~data.badChannels, 'nans');
tsBlankAverage = to157chan(tsBlankAverage', ~data.badChannels, 'nans');

tsFullBlankDiff = tsFullAverage-tsBlankAverage;

% Define time points to plot
onecyclems = (1/12) * 1000;
nrCyclesToPlot = 5;
timePoints = 1:1:ceil(nrCyclesToPlot*onecyclems);

%% Visualize predictions
% Define image size
scr_sz = get(0,'ScreenSize');
rows = ceil(scr_sz(3)/3);
cols = ceil(scr_sz(3)/3);

fH = figure(1); clf; set(gcf, 'Position', [1 1 rows, cols], 'Color', 'w');

% Set colormap, colormap limits
cmap     = bipolar;
clim     = [-1,1].*abs(min(tsFullBlankDiff(:)));
fprintf('\n(%s): Creating frames.', mfilename);

F(length(timePoints)) = struct('cdata',[], 'colormap', []);
for t = 1:length(timePoints)
    cla;
    [~,ch] = megPlotMap(tsFullBlankDiff(timePoints(t),:),clim,[],cmap,sprintf('t = %3.0f ms',timePoints(t)),[],[], ...
        'mask', 'convex', ...
        'pointsize', 10);
    set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12);
    drawnow;
    
    F(t) = getframe(fH);
end
    
 
% Create movie
mov = implay(F, 4);

if saveFig
    figureDir = fullfile(fmsRootPath,'figures',subject{subjectsToPlot});
    fname = fullfile(figureDir, 'singleCycleMovieDataFigure3');
    vid = VideoWriter(fname, 'MPEG-4');
    open(vid);
    
    writeVideo(vid, F);
end

fprintf('\n(%s): Done!', mfilename)


