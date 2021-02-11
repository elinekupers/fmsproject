%% Script to plot Figure 5, Data
%
% Winawer, Kay, Foster, Parvizi, Wandell.
% *Asynchronous broadband signals are the principal source of the BOLD
% response in human visual cortex*
% _Current Biology, 2013_
%
% This script plots the times series and the spectal responses from an ECoG
% channel to large-field flicker. It produces the two panels on the left
% side of figure 5. Data comes from subject 1, channel 66 (V1/V2 fovea).
%
%
% Code by Jonathan Winawer, 2013
%
% (c) 2013, Stanford vistalab
%

%% Set some parameters

savepth   = fullfile(ecogPRFrootPath, 'scratch');
calcPower = true;  % Plot spectra as power (squared amplitude) 
useHann   = false; % Use a Hanning window for spectral analysis 
fmax      = 150; 
stimF     = 15;

%% Load the data

% This includes  
%   t:          time vector (seconds), 1x3 cell for 3 runs
%   ts:         raw time series (microvolts), 1x3 cell for 3 runs
%   onsets:     epoch onsets in temporal samples (not seconds), 1x3 cell
%   sampleRate: ECoG sampling rate, in Hz
%   T:          epoch length (in seconds)
%   subjnum:    subject number (corresponds with numbering in paper)
%   runType:    indicates that this data comes from OnOff expts
%   dataType:   indicates that data was referenced to common average
%
%   Note that each run consisted of 6 'on' epochs, followed by 6
%   'off' epochs, repeated 4 times (i.e., 4 on-off blocks of duration 12*T
%   seconds each)

load(fullfile(ecogPRFrootPath, 'data', 'figure5Data'));

%% Compute spectra of each epoch

% We have 3 on-off runs
numruns = numel(ts);

tmp = cell(1,numruns);
for run = 1:numruns
    tmp{run} = ecogTSeriesVector2TSeriesMatrix(ts{run}, onsets{run});
end

tsmatrix = catcell(2, tmp); clear tmp

epochs    = reshape(1:size(tsmatrix,2), 6, []);
onepochs  = epochs(:,1:2:end); onepochs = onepochs(:);
offepochs = epochs(:,2:2:end); offepochs = offepochs(:);

on.signal  = tsmatrix(:,onepochs);
off.signal = tsmatrix(:,offepochs);
%% Summarize spectra and time series
epochT = linspace(0, T, size(tsmatrix,1)+1); epochT = epochT(1:end-1);

% compute means across trials
[on, off] = ecogCalcOnOffSpectra(on, off, useHann, calcPower);

% Plot time series and spectra
fH = ecogPlotOnOffSpectra(on, off, epochT, stimF, calcPower);

%% Save

hgexport(fH(1), fullfile(savepth, 'Figure5_Data_TimeSeries.eps'));
hgexport(fH(2), fullfile(savepth, 'Figure5_Data_Spectra.eps'));

%%
% % Note for author
% % This is how the saved data file was made.
% 
% chan = 66; 
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
% runType = 'onOff';
% dataType = 'Common Average Rereferenced';
% save ~/matlab/git/ECoG_PURL/data/figure5Data.mat t ts onsets T sampleRate subjnum runType dataType
% 
% 

