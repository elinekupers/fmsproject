function plotSpectrumIndividualChannel(exampleSubject)

% Function to plot the spectrum of an individual channel for an individual
% subject.

subject         = {'wlsubj002','wlsubj004','wlsubj005','wlsubj006','wlsubj010','wlsubj011','wlsubj048', 'wlsubj046','wlsubj039','wlsubj059', 'wlsubj067', 'wlsubj070'};
if nargin < 1; exampleSubject  = 12; end % Which example subject to show if not defined

% Set up paths
figureDir              = fullfile(fmsRootPath, 'figures', subject{exampleSubject}); % Where to save images?
dataDir                = fullfile(fmsRootPath, 'data');    % Where to get data?
saveFigures            = true;      % Save figures in the figure folder?


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
colors    = [0 147 68; 126 126 126]./255;

% Spectrum of example channel
fH = figure('position',[0,300,500,500]); clf(fH); set(fH, 'Name', 'Spectrum of one MEG sensor' , 'NumberTitle', 'off');

% Define axes
f = (0:999);
xl = [8 150];
fok = f;
fok(f<=xl(1) | f>=xl(2) | mod(f,60) < 2 | mod(f,60) > 58 ) = [];
xt = [12:12:72, 96,144];
yt = -29:-25;
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
plot(fok, mn.full(:,2),  '-',  'Color', colors(1,:), 'LineWidth', 2); hold on;
plot(fok, mn.blank(:,2),  '-',  'Color', colors(2,:), 'LineWidth', 2);

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
plot(fok, mn.full(:,2),  '-',  'Color', colors(1,:), 'LineWidth', 2); hold on;
plot(fok, mn.blank(:,2),  '-',  'Color', colors(2,:), 'LineWidth', 2);

% Reset x scale and y tickmarks
xl = [60 150];
yt = [-29 -27.8];

% format x and y axes
set(gca, 'XLim', xl, 'XTick', xt, 'XScale', 'log', 'YScale','log');
set(gca,'ytick',10.^yt, 'ylim',10.^yt, 'TickDir', 'out', 'FontSize', 12);
box off;

% label figure, add stimulus harmonic lines, and make it look nice
xlabel('Frequency (Hz)'); ylabel('Power (T^2)');
title(sprintf('Channel %d', sensorIdx));
yl2 = get(gca, 'YLim');
for ii =12:12:180, plot([ii ii], yl2, '-', 'Color', [100 100 100]./255); end

if saveFigures
    figurewrite(fullfile(figureDir, 'Figure2A_SpectrumZoomBroadband'), [],0,'.',1);
end


%% Visualize where channel is on mesh
figure;
sensorloc = ones(1,size(data.ts,1))*0.5;
sensorloc(sensorIdx)=1;
megPlotMap(to157chan(sensorloc,~data.badChannels,'zeros'),[0 1],[],'gray', [], [], [], 'interpmethod', 'nearest'); colorbar off


if saveFigures
    figurewrite(fullfile(figureDir, 'Figure2A_LocationOfExampleChannel'), [],0,'.',1);
end


return
