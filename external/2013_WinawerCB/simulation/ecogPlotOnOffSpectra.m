function fH = ecogPlotOnOffSpectra(on, off, t, stimF, calcPower)
% Plot time series and spectra from on-off experiments
% 
%   fH = ecogPlotOnOffSpectra(on, off, t, stimF, calcPower)
%
% Inputs
%   on: struct containing a field called 'signal', which is a matrix 
%       num time points x epochs
%   off: same, but for off epochs
%   t:  a vector of sample times for one epoch (in seconds)
%   stimF: frequency of stimulus locked responses
%   calcPower: (boolean) If true, plot spectra as squared amplitude
%
% Outputs
%   fH:  vector of figure handles. 
%   fH(1) is handle to spectral plot
%   fH(2) is handle to time series plot
%
%
% Example
%   t = .001:.001:1;
%   ntrials = 50;
%   on.signal  = 75*(randn(length(t),ntrials) + sin(repmat(t',1, ntrials)*2*pi*10));
%   off.signal = 75*(randn(length(t),ntrials)) ;
%   [on off] = ecogCalcOnOffSpectra(on, off);
%   fH = ecogPlotOnOffSpectra(on, off, t, 10, true)

fmax = 150;
if notDefined('calcPower'), calcPower = false; end
if size(t, 1) == 1, t = t'; end
f = (0:length(t)-1)/max(t);
f = f(f<=fmax);

%% PLOT SPECTRUM
fH(1) = figure;

set(gcf, 'Color', 'w')

cols.plot = [0 0 0; 126 126 126]./255;

xticks = stimF:stimF:150;
xtickLabels = cellstr(num2str(ceil(xticks')))';
xtickLabels{5} = '';
xtickLabels{7} = '';
xtickLabels{8} = '';


if calcPower, yl = [10^-2 10^2].^2; else yl = [10^-2 10^2]; end
set(gca, 'Xlim', [8 fmax], 'XScale', 'log', 'YScale', 'log',...
    'XTick', xticks,  'XTickLabel', xtickLabels,  ...
    'YLim', yl, 'FontSize', 16, ...
    'ColorOrder', cols.plot)
hold all


fs = 1:length(f);
plot(...
    f, on.meanFFT(fs), '-', f, off.meanFFT(fs), '-', 'LineWidth',2)...
%     f, on.fftMean(fs), '-', f, off.fftMean(fs), '-', ...
%     'LineWidth', 2)

set(gca, 'XGrid', 'on', 'GridLineStyle', '-', 'XMinorGrid', 'off');
legend({'Mean spectrum ON' 'Mean spectrum OFF'}, ...
    'Location', 'Southwest',  'Box', 'off', 'EdgeColor', 'w', 'Color', 'none')
%     'Spectrum of mean time series ON' 'Spectrum of mean time series OFF'},...
    

xlabel('Frequency (Hz)')
if calcPower, ylabel('Power (µV^2)'); else ylabel('Amplitude (µV)'); end


%% PLOT TIME SERIES
fH(2) = figure;
set(gcf, 'Color', 'w')

xt = (0:stimF*max(t))/stimF;
xtLabels = cellstr(num2str(round(xt'*100)/100))';
xtLabels(1:2:end) = {' '};


set(gca, 'ColorOrder', cols.plot, 'XTick', xt, ...
    'XTickLabel', xtLabels, ...
    'XLim', [0 max(t)],  'XGrid', 'on', ...
    'FontSize', 16, 'YLim', [-100 100], 'YTick', -100:50:100);
hold all;


plot(t, on.meanTs, t, off.meanTs, 'LineWidth', 2);
xlabel('Time (s)')
ylabel('Amplitude (µV)')


return

