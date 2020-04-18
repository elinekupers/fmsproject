function makeFigure2(varargin)
% Function to plot the spectrum of an individual channel for an individual
% subject for the manuscript:
%   A visual encoding model links magnetoencephalography
%   signals to neural synchrony in human cortex.
%       by Kupers, Benson, Winawer (YEAR) JOURNAL.
%
% INPUTS:
%   [subjectsToPlot]        :  (int)  subject nr to plot (default: 12)
%
% Example: Plot example subject in manuscript (S12)
%  makeFigure2('subjectsToPlot', 12, 'saveFig', true)
%
%
% By Eline Kupers, NYU (2017)
p = inputParser;
p.KeepUnmatched = true;
p.addParameter('subjectsToPlot', 12);
p.addParameter('saveFig', true, @islogical);
p.parse(varargin{:});

% Rename variables
subjectsToPlot        = p.Results.subjectsToPlot;
saveFig               = p.Results.saveFig;

% Get subject names and their corresponding data session number
[subject, dataSession] = getSubjectIDs;

% Set up paths
figureDir              = fullfile(fmsRootPath, 'figures', subject{exampleSubject}); % Where to save images?
dataDir                = fullfile(fmsRootPath, 'data');    % Where to get data?


%% 1. Load subject's data

% Go from subject to session nr
whichSession = dataSession(subjectsToPlot);
    
% Get amplitude data
data = loadData(fullfile(dataDir, subject{exampleSubject}), whichSession,'timeseries');

% Define MEG sensor to plot
sensorIdx = 1;
               
% Define color to plot full conditions
colors    = [0 0 0; 126 126 126]./255;

% Define axes
f = (0:999);
xl = [8 150];
fok = f;
fok(f<=xl(1) | f>=xl(2) | mod(f,60) < 2 | mod(f,60) > 58 ) = [];
xt = [12:12:72, 96,144];

% compute spectrum
spec = abs(fft(squeeze(data.ts(sensorIdx, :,:))))/size(data.ts,2)*2;

% compute power for epochs corresponding to a condition and
% trim data to specific frequencies
dataFull = spec(:,data.condEpochsFull).^2;
dataFull = dataFull(fok+1,:);

dataBlank = spec(:,data.condEpochsBlank).^2;
dataBlank = dataBlank(fok+1,:);

% Get the median
mn.full = prctile(dataFull,[16,50,84],2);
mn.blank = prctile(dataBlank,[16,50,84],2);

% Multiply with 10^15 to get values in units of fT
if any(intersect(whichSession, 9:14))
    mn.full = mn.full.*10^-15.*10^-15;
    mn.blank = mn.blank.*10^-15.*10^-15;
end

%% Spectrum of example channel
fH = figure('position',[0,300,500,500]); clf(fH); set(fH, 'Name', 'Spectrum of one MEG sensor' , 'NumberTitle', 'off');

% plot median 
plot(fok, mn.full(:,2),  '-',  'Color', colors(1,:), 'LineWidth', 2); hold on;
plot(fok, mn.blank(:,2),  '-',  'Color', colors(2,:), 'LineWidth', 2);

% format x and y axes
yt = [1:5];
yl=[yt(1),yt(end)];
set(gca, 'XLim', xl, 'XTick', xt, 'XScale', 'log', 'YScale','log');
set(gca,'ytick',10.^yt, 'ylim',10.^yl, 'TickDir', 'out', 'FontSize', 12);
box off;

% label figure, add stimulus harmonic lines, and make it look nice
xlabel('Frequency (Hz)'); ylabel('Power (fT^2)');
title(sprintf('Channel %d', sensorIdx));
yl2 = get(gca, 'YLim');
for ii =12:12:180, plot([ii ii], yl2, '-', 'Color', [100 100 100]./255); end

if saveFig
    figurewrite(fullfile(figureDir, 'Figure2A_SpectrumOneChannel'), [],0,'.',1);
end

%% Plot zoom into broadband frequencies
fH2 = figure('position',[567, 655, 300, 281]); clf(fH2); set(fH2, 'Name', 'Spectrum of one MEG sensor' , 'NumberTitle', 'off');

% plot median
plot(fok, mn.full(:,2),  '-',  'Color', colors(1,:), 'LineWidth', 2); hold on;
plot(fok, mn.blank(:,2),  '-',  'Color', colors(2,:), 'LineWidth', 2);

% Reset x scale and y tickmarks
xl = [60 150];
yl = [1 2.1];
yt = [1 2];

% format x and y axes
set(gca, 'XLim', xl, 'XTick', xt, 'XScale', 'log', 'YScale','log');
set(gca,'ytick',10.^yt, 'ylim',10.^yl, 'TickDir', 'out', 'FontSize', 12);
box off;

% label figure, add stimulus harmonic lines, and make it look nice
xlabel('Frequency (Hz)'); ylabel('Power (T^2)');
title(sprintf('Channel %d', sensorIdx));
yl2 = get(gca, 'YLim');
for ii =12:12:180, plot([ii ii], yl2, '-', 'Color', [100 100 100]./255); end

if saveFig
    figurewrite(fullfile(figureDir, 'Figure2A_SpectrumZoomBroadband'), [],0,'.',1);
end


%% Visualize where channel is on mesh
figure;
sensorloc = ones(1,size(data.ts,1))*0.5;
sensorloc(sensorIdx)=1;
megPlotMap(to157chan(sensorloc,~data.badChannels,'zeros'),[0 1],[],'gray', [], [], [], 'interpmethod', 'nearest'); colorbar off
title('White = Sensor location, Black = bad sensor');

if saveFigures
    figurewrite(fullfile(figureDir, 'Figure2Inset_LocationOfExampleChannel'), [],0,'.',1);
end


%% Calculate increase in SL and BB

% Define frequencies to compute the broadband power
fs           = 1000;         % Sample rate
f            = 0:150;        % Limit frequencies to [0 150] Hz
slFreq       = 12;           % Stimulus-locked frequency
tol          = 1.5;          % Exclude frequencies within +/- tol of sl_freq
slDrop       = f(mod(f, slFreq) <= tol | mod(f, slFreq) > slFreq - tol);

% Exclude all frequencies below 60 Hz when computing broadband power
lfDrop       = f(f<60);

% Define the frequenies and indices into the frequencies used to compute
% broadband power
[~, abIndex] = setdiff(f, [slDrop lfDrop]);

% Create function handles for the frequencies that we use
keepFrequencies    = @(x) x(abIndex);

for iter = 1:1000
    fullIdx = find(data.condEpochsFull==1);
    subsetFullEpochs = fullIdx(randi([1, length(fullIdx)],1, sum(data.condEpochsBlank)));

    % Get data
    dataIn = cat(1, data.ts(sensorIdx,:,subsetFullEpochs), data.ts(sensorIdx,:,data.condEpochsBlank));
    bbFull  = log10(getbroadband(dataIn(1,:,:),keepFrequencies,fs)); % Broadband function gives data in units of power
    bbBlank = log10(getbroadband(dataIn(2,:,:),keepFrequencies,fs)); % Broadband function gives data in units of power
        
    slFull  = log10(getstimlocked(dataIn(1,:,:), 13).^2);  % Square to get units of power
    slBlank = log10(getstimlocked(dataIn(2,:,:), 13).^2); % Square to get units of power
              
    % Full field minus blank, averaged across top 5 channels
    diffBB = mean(bbFull - bbBlank);
        
    % Add individual subjects stimulus locked results to all results
    diffSL = mean(slFull - slBlank);

    % Calculate percent change
    percentdiff.bb(iter,:) = 100*((10.^diffBB)-1);
    percentdiff.sl(iter,:) = 100*((10.^diffSL)-1);
end

fprintf('Mean change in broadband power: %4.1f%%, with sd %1.1f%%\n', mean(percentdiff.bb),std(percentdiff.bb));
fprintf('Mean change in stimulus locked power: %4.1f%%, with sd %1.1f%%\n', mean(percentdiff.sl),std(percentdiff.sl));

return
