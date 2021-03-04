function makeFigure3PhaseMaps(varargin)
% This is a function to make phase maps, instead of amplitude maps,
% similar to Figure 3 from the manuscript:
%   A visual encoding model links magnetoencephalography
%   signals to neural synchrony in human cortex.
%       by Kupers, Benson, Winawer (YEAR) JOURNAL.
%
% This figure shows the phase maps for the stimulus-locked response and 
% an asynchronous broadband response to a large field flickering (12 Hz) 
% dartboard pattern.
%
% To runs this script, you need:
% (1) the data from the denoiseproject in the data folder of its FMS
%     code repository
%
% (2) MEG_utils toolbox added to the paths. For example:
%     tbUse('ForwardModelSynchrony');
%        or to only add the MEG_utils toolbox:
%     addpath(genpath('~/matlab/git/toolboxes/meg_utils'))
%
% INPUTS:
%   [subjectsToPlot]     : (int)  subject nr to plot (default: 12)
%   [plotMeanSubject]    : (bool) plot average across all 12 subjets or not?
%                                 (default: true)
%                                 
%
% Example 1: Plot first subject
%  makeFigure3PhaseMaps('subjectsToPlot', 1, 'plotMeanSubject', false, 'saveFig', true)
% Example 2: Plot example subject in manuscript (S12)
%  makeFigure3PhaseMaps('subjectsToPlot', 12, 'plotMeanSubject', false, 'saveFig', true)
% Example 3: Plot all subjects and group average
%  makeFigure3PhaseMaps('subjectsToPlot', 1:12, 'plotMeanSubject', true, 'saveFig', true)
%
% By Eline Kupers (NYU) 2021

%% 0. Set up paths and define parameters
p = inputParser;
p.KeepUnmatched = true;
p.addParameter('subjectsToPlot', 12);
p.addParameter('plotMeanSubject', true, @islogical); 
p.parse(varargin{:});

% Rename variables
subjectsToPlot          = p.Results.subjectsToPlot;
plotMeanSubject         = p.Results.plotMeanSubject;

% Get subject names and corresponding data session number
[subject, dataSession] = getSubjectIDs;

% Set up paths
figureDir        = fullfile(fmsRootPath, 'figures'); % Where to save images?
dataDir          = fullfile(fmsRootPath, 'data');    % Where to get data?

% Plotting params
markerType    = '.';

% Load all subjects when plotting the mean
if plotMeanSubject
    subjectsToLoad = 1:length(subject);
else
    subjectsToLoad = subjectsToPlot;
end

allData = cell(size(subjectsToLoad));

slIdx = 12+1;

%% 1. Load subject's data
for s = subjectsToLoad
    
    % Go from subject to session nr
    whichSession = dataSession(s);
    
    % Load data
    data     = loadData(fullfile(dataDir, subject{s}), whichSession, 'type', 'timeseries');
    dataFull = data.ts(:,:,data.condEpochsFull==1);
    
    % Average time series across full field epochs
    tsFullAverage = mean(dataFull,3);

    % Fourier transform Full field data, take phase at 12 Hz
    spec = fft(tsFullAverage,[],2);
    slAmpsFull   = abs(squeeze(spec(:,slIdx,:)))' / size(data,2)*2;
    slPhaseFull  = angle(squeeze(spec(:,slIdx,:)))' / size(data,2)*2;
   
    % Add bad channels back in as nans
    slAmpsFull  = to157chan(slAmpsFull, ~data.badChannels, 'nans');
    slPhaseFull = to157chan(slPhaseFull, ~data.badChannels, 'nans');
    
    allData{s}.slAmpsFull = slAmpsFull;
    allData{s}.slPhaseFull = slPhaseFull;
    
end

%% 3. Mask data by SNR, plot data per subject

if length(subjectsToPlot) >1
    rows = 2;
    cols = ceil(length(subjectsToPlot)/2);
else
    rows = 1;
    cols = 1;
end

figure(1); clf; set(gcf, 'color','w','Name','Observed 12 Hz Amplitude Data')
figure(2); clf; set(gcf, 'color','w','Name','Observed 12 Hz Phase Data')

for s = subjectsToPlot

    unreliableSensors = allData{s}.slAmpsFull > prctile(allData{s}.slAmpsFull,80);
    allData{s}.slAmpsFull(~unreliableSensors) = NaN;
    allData{s}.slPhaseFull(~unreliableSensors) = NaN;
    
    dataToPlot  = cat(1, allData{s}.slAmpsFull, ...
        allData{s}.slPhaseFull);

    figure(1);
    subplot(rows,cols,s);  hold all;
    megPlotMap(dataToPlot(1,:), [],[],[],[],[],[],'interpmethod','nearest'); colorbar;
    title(sprintf('S%d', s)) 
    
    figure(2);
    subplot(rows,cols,s); hold all; 
    megPlotMap(dataToPlot(2,:), [-pi pi],[],'hsv',[],[],[],'interpmethod','nearest'); colorbar;
    title(sprintf('S%d', s)) 
    
end


%% 4. Plot average across subjects if requested

if plotMeanSubject
    
    for ii = subjectsToLoad
        phase(ii,:) = allData{ii}.slPhaseFull;
        amps(ii,:) = allData{ii}.slAmpsFull;
    end
    
    mnPhase = NaN(1,size(phase,2));
    for sensor = 1:size(phase,2)
        nanIdx = isnan(phase(:,sensor));
        mnPhase(1,sensor) = circ_mean(phase(~nanIdx,sensor),[],1);
    end
    
    mnAmps = mean(amps,1,'omitnan');
    
    unreliableSensors = mnAmps > prctile(mnAmps,80);
    mnAmps(~unreliableSensors) = NaN;
    mnPhase(~unreliableSensors) = NaN;
    dataToPlot  = cat(1, mnAmps,mnPhase);

    figure(3); clf
    megPlotMap(dataToPlot(1,:), [],[],[],[],[],[],'interpmethod','nearest'); colorbar;
    title('12 Hz Amplitude Group Average (N=12)') 
    
    figure(4); clf
    megPlotMap(dataToPlot(2,:), [-pi pi],[],'hsv',[],[],[],'interpmethod','nearest'); colorbar;
    title('12 Hz Phase Group Average (N=12)') 

end

return



