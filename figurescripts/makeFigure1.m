function makeFigure1(exampleSubject)

% Function to plot the spectrum of an individual channel for an individual
% subject.

subject         = {'wlsubj002','wlsubj004','wlsubj005','wlsubj006','wlsubj010','wlsubj011','wlsubj048', 'wlsubj046','wlsubj039','wlsubj059', 'wlsubj067', 'wlsubj070'};
if nargin < 1; exampleSubject  = 12; end % Which example subject to show if not defined

% Set up paths
figureDir              = fullfile(fmsRootPath, 'figures', subject{exampleSubject}); % Where to save images?
dataDir                = fullfile(fmsRootPath, 'data');    % Where to get data?
saveFigures            = false;      % Save figures in the figure folder?


%% 1. Load subject's data

switch subject{exampleSubject}
    % Go from subject to session nr
    case 'wlsubj002'
        whichSession = 2;
    case 'wlsubj004'
        whichSession = 7;
    case 'wlsubj005'
        whichSession = 8;
    case 'wlsubj006'
        whichSession = 1;
    case 'wlsubj010'
        whichSession = 6;
    case 'wlsubj011'
        whichSession = 5;
    case 'wlsubj048'
        whichSession = 9; % Full field Only
    case 'wlsubj046'
        whichSession = 10; % Full field Only
    case 'wlsubj039'
        whichSession = 11; % Full field Only
    case 'wlsubj059'
        whichSession = 12; % Full field Only
    case 'wlsubj067'
        whichSession = 13; % Full field Only
    case 'wlsubj070'
        whichSession = 14; % Full field Only
end
    
% Get amplitude data
data = loadData(fullfile(dataDir, subject{exampleSubject}), whichSession,'timeseries');

% Define MEG sensor to plot
sensorIdx = 1;
               
% Define color to plot full conditions
colors    = [0 0 0; 126 126 126]./255;

% Spectrum of example channel
fH = figure('position',[0,300,500,500]); clf(fH); set(fH, 'Name', 'Spectrum of one MEG sensor' , 'NumberTitle', 'off');

% Define axes
f = (0:999);
xl = [8 150];
fok = f;
fok(f<=xl(1) | f>=xl(2) | mod(f,60) < 2 | mod(f,60) > 58 ) = [];
xt = [12:12:72, 96,144];
yt = 1:5;
yl=[yt(1),yt(end)];

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

% plot median
plot(fok, mn.full(:,2).*10^15.*10^15,  '-',  'Color', colors(1,:), 'LineWidth', 2); hold on;
plot(fok, mn.blank(:,2).*10^15.*10^15,  '-',  'Color', colors(2,:), 'LineWidth', 2);

% format x and y axes
set(gca, 'XLim', xl, 'XTick', xt, 'XScale', 'log', 'YScale','log');
set(gca,'ytick',10.^yt, 'ylim',10.^yl, 'TickDir', 'out', 'FontSize', 12);
box off;

% label figure, add stimulus harmonic lines, and make it look nice
xlabel('Frequency (Hz)'); ylabel('Power (fT^2)');
title(sprintf('Channel %d', sensorIdx));
yl2 = get(gca, 'YLim');
for ii =12:12:180, plot([ii ii], yl2, '-', 'Color', [100 100 100]./255); end

if saveFigures
    figurewrite(fullfile(figureDir, 'Figure2A_SpectrumOneChannel'), [],0,'.',1);
end

%% Plot zoom into broadband frequencies
fH2 = figure('position',[567, 655, 300, 281]); clf(fH2); set(fH2, 'Name', 'Spectrum of one MEG sensor' , 'NumberTitle', 'off');

% plot median
plot(fok, mn.full(:,2).*10^15.*10^15,  '-',  'Color', colors(1,:), 'LineWidth', 2); hold on;
plot(fok, mn.blank(:,2).*10^15.*10^15,  '-',  'Color', colors(2,:), 'LineWidth', 2);

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

if saveFigures
    figurewrite(fullfile(figureDir, 'Figure1A_SpectrumZoomBroadband'), [],0,'.',1);
end


%% Visualize where channel is on mesh
figure;
sensorloc = ones(1,size(data.ts,1))*0.5;
sensorloc(sensorIdx)=1;
megPlotMap(to157chan(sensorloc,~data.badChannels,'zeros'),[0 1],[],'gray', [], [], [], 'interpmethod', 'nearest'); colorbar off


if saveFigures
    figurewrite(fullfile(figureDir, 'Figure1A_LocationOfExampleChannel'), [],0,'.',1);
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
percentdiff.bb = 100*((10.^diffBB)-1);
percentdiff.sl = 100*((10.^diffSL)-1);

fprintf('Mean change in broadband power: %4.1f%% \n', percentdiff.bb);
fprintf('Mean change in stimulus locked power: %4.1f%% \n', percentdiff.sl);

return
