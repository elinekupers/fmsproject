%% Script to plot Figure 1c  
% 
% Winawer, Kay, Foster, Parvizi, Wandell.
% *Asynchronous broadband signals are the principal source of the BOLD
% response in human visual cortex*
%  _Current Biology, 2013_
%
% This figure shows an example time series from an On and Off flickering
% large-field contrast pattern. The flicker was 7.5 Hz square wave
% (contrast reversals 15 times per second). The subject was S1 and the
% channel 104 (V1/V2 periphery). 
%
% Copyright Jonathan Winawer, 2013

%% Set up paths and parameters


savepth  = fullfile(ecogPRFrootPath, 'scratch');
datafile = fullfile(ecogPRFrootPath, 'data', 'figure1Data');
run      = 2; % There are several on-off time series for this subject. 
              % We plot the time series from run 2.


%% Load the data

% This includes  
%   t:          time vector (seconds), 1x3 cell for 3 runs
%   ts:         raw time series (microvolts), 1x3 cell for 3 runs
%   onsets:     epoch onsets in temporal samples, 1x3 cell
%   sampleRate: ECoG sampling rate, in Hz
%   T:          epoch length (in seconds)
%   subjnum:    subject number (corresponds with numbering in paper)
%   runType:    indicates that this data comes from OnOff expts
%   dataType:   indicates that data was re-referenced to common average
%
%   Note that the experiment consisted of 6 'on' epochs, followed by 6
%   'off' epochs, repeated 4 times (i.e., 4 on-off blocks of duration 12*T
%   seconds each)

load(datafile);

% We have 3 runs of the same type. We will plot one time series (run 2).
t       = t{run};
ts      = ts{run};
onsets  = onsets{run};

% The last sample is 1 epoch length after the last epoch onset
lastsample  = onsets(end) + round(T*sampleRate)-1;
firstsample = onsets(1);

% For purposes of plotting, we will show 6 consecutive ON epochs, followed
% by 6 consecutive OFF epochs, repeated 4 times. 
onsets =  onsets(1:6:end);
offsets = [onsets(2:end) lastsample];

%% Plot the whole time series

% Set up the figure 
fH = figure; clf; set(fH, 'Color', 'w'); hold on

yl = std(ts) * [-5 5];
xl = [0 t(lastsample)];

set(gca, 'XLim', xl, 'XTick', 0:12:48, ...
    'YLim', yl, 'YTick', -200:100:200, 'FontSize', 16)

xlabel('Time (s)')
ylabel('Voltage (µV)')

% Plot the entire time series in light gray
plot(t(firstsample:lastsample), ts(firstsample:lastsample),  'Color', [.7 .7 .7]);

%% Plot the OFF epochs in dark gray
figure(fH)

% loop through the OFF epochs 
for ii = 2:2:8
    thesesamples = onsets(ii):offsets(ii);
    plot(t(thesesamples), ts(thesesamples),  'Color', [.3 .3 .3]);
end

hgexport(fH,  fullfile(savepth, 'Figure1C_onOffTS.eps'));


return



%% 
% % Note for author
% % This is how the saved data file was made.
% 
% chan = 104; 
% s = 9; % note that s9 is subj 1 in the manuscript
% for run = 1:3; 
%     t{run}  = ecogGet(s, 't', 'onoff', run, chan); 
%     ts{run} = ecogGet(s, 'ts', 'onoff', run, chan, 'CAR'); 
%     onsets{run}  = ecogGet(s, 'trial onsets', 'onoff', run, 'frames');
% end
% 
% T = ecogGet(s, 'tr');
% sampleRate = ecogGet(s, 'fs');
% 
% subjnum = 1;
% channel = 104;
% runType = 'onOff';
% dataType = 'Common Average Rereferenced';
% save ~/matlab/git/ECoG_PURL/data/figure1Data.mat t ts onsets T sampleRate subjnum runType dataType






