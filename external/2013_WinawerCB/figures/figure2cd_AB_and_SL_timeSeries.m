%% Script to plot Figure 2c,d  
%
% Winawer, Kay, Foster, Parvizi, Wandell
% *Asynchronous broadband signals are the principal source of the BOLD
% response in human visual cortex*
% _Current Biology, 2013_
%
% This figure shows the Asynchronous Broadband (panel c) and
% Stimulus-Locked (panel d) time series in response to wide and narrow bar
% stimuli. The data is from subject 2, channel 65 (V1/V2 periphery).
%
% Copyright Jonathan Winawer, 2013


%% Set up paths and parameters

savepth   = fullfile(ecogPRFrootPath, 'scratch');
datafile  = fullfile(ecogPRFrootPath, 'data', 'figure2Data');
fmax      = 150;  % spectral maximum (Hz) for fittiing broadband time series
useHann   = true; % apply a Hanning window before computing spectrum

%% Load the data

% This includes  
%   t:          time vector (seconds), 1x3 cell for 3 runs
%   ts:         raw time series (microvolts), 1x3 cell for 3 runs
%   onsets:     sample numbers indicating epoch onsets, 1x3 cell
%   sampleRate: ECoG sampling rate, in Hz
%   T:          epoch length (in seconds)
%   subjnum:    subject number (corresponds with numbering in paper)
%   runType:    indicates the type of stimuli
%   dataType:   indicates that data was referenced to common average

load(datafile);

%% Compute asynchronous broadband (ab) and stimulus-locked (sl) time series

% Convert the time series to spectra (384 epochs x 151 frequencies)
spectra = ecogGetSpectra(ts, onsets, T, fmax, useHann);

% Summarize spectra as broadband and stimulus-locked timeseries (384x1 each)
spectralSummary = ecogSummarizeSpectra(spectra, [], T, fmax);

%  Reshape the time series into a matrix of time points (96) by runs (4)
ab.all = reshape(spectralSummary.ab, 96,4); 
sl.all = reshape(spectralSummary.sl, 96,4);

% Average across runs of the same type
ab.average = [mean(ab.all(:,1:2),2) mean(ab.all(:,3:4),2)];
sl.average = [mean(sl.all(:,1:2),2) mean(sl.all(:,3:4),2)];

% normalize so that the mean response during stimulus blanks is zero
ab.normalized = ecogSubtractBlanks(ab.average, runType);
sl.normalized = ecogSubtractBlanks(sl.average, runType);

%% Plot
fH = figure; pos = get(fH, 'pos'); 
set(fH, 'Color', 'w', 'position', [pos(1) pos(2) 600 800]);

% Plot the broadband time series
subplot(2,1,1)
set(gca, 'ColorOrder', [0 0 0; 1 0 0], 'FontSize', 16, ...
    'XTick', 0:12:96, 'XLim', [0 96], 'XGrid', 'on'); hold on
xlabel('Time (s)'); ylabel('Power (µV^2)')

plot(1:96, ab.normalized, 'LineWidth', 2)

title('Broadband time series')
legend('Narrow bars', 'Wide bars')

subplot(2,1,2)
set(gca, 'ColorOrder', [0 0 0; 1 0 0], 'FontSize', 16, ...
    'XTick', 0:12:96, 'XLim', [0 96], 'XGrid', 'on'); hold on
xlabel('Time (s)'); ylabel('Amplitude (µV)')

plot(1:96, sl.normalized, 'LineWidth', 2)

title('Stimulus-locked time series')
legend('Narrow bars', 'Wide bars')


%% SAVE

hgexport(fH, fullfile(savepth, 'Figure2cd_AB_and_SL_timeSeries.eps'));

return



%% Note for author
% %
% % This is how the saved data file was made.
%
% chan = 65;
% s = 13; % note that s13 is subj 2 in the manuscript
% runType = {'barThin', 'barThick'};
% for ii = 1:length(runType)
%     exp = runType{ii};
%     for run = ecogGet(s, 'runs', exp);
%         t.(exp){run}  = ecogGet(s, 't', exp, run, chan);
%         ts.(exp){run} = ecogGet(s, 'ts', exp, run, chan, 'CAR');
%         onsets.(exp){run}  = ecogGet(s, 'trial onsets', exp, run, 'frames');
%     end
% end
% T = ecogGet(s, 'tr');
% sampleRate = ecogGet(s, 'fs');
% 
% subjnum = 2;
% dataType = 'Common Average Rereferenced';
% save ~/matlab/git/ECoG_PURL/data/figure2Data.mat t ts onsets ...
%     T sampleRate subjnum runType dataType
