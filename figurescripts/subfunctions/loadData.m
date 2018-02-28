function [data, badChannels] = loadData(dataDir, whichSubject)

%% Define parameters to get SL and BB data
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



%% Load data
load(sprintf(fullfile(dataDir, 's%02d_conditions.mat'),whichSubject));
load(sprintf(fullfile(dataDir, 's%02d_sensorData.mat'),whichSubject));
load(sprintf(fullfile(dataDir, 's%02d_denoisedData_bb.mat'),whichSubject));
load(sprintf(fullfile(dataDir, 's%02d_denoisedts.mat'),whichSubject));


% preprocessing parameters (see nppPreprocessData)
varThreshold        = [0.05 20];
badChannelThreshold = 0.2;
badEpochThreshold   = 0.2;
dataChannels        = 1:157;

% Preprocess raw sensordata
[sensorData, badChannels0, badEpochs0] = nppPreprocessData(sensorData(:,:,dataChannels), ...
    varThreshold, badChannelThreshold, badEpochThreshold, false);

% ---- Define first epochs in order to remove later ------------------
badEpochs0(1:6:end) = 1;

% Remove bad channels and bad epochs from data and conditions
sensorData = sensorData(:,~badEpochs0, ~badChannels0);

% Permute sensorData for denoising
sensorData = permute(sensorData, [3 1 2]);

% time domain data before and after denoising
denoisedts = denoisedts_bb;

sl_ts = sensorData;
bb_ts = denoisedts{1};

design = zeros(size(conditions,1),3);
design(conditions == 1,1) = 1; % Full
design(conditions == 5,2) = 1; % Right
design(conditions == 7,3) = 1; % Left

% Get rid of bad epochs
design_sl = design(~badEpochs0,:);
design_bb = design(~badEpochs,:);

% Define conditions: Full, right, left, off
condEpochsSL = {design_sl(:,1)==1, design_sl(:,2)==1, design_sl(:,3)==1, all(design_sl==0,2)};
condEpochsBB = {design_bb(:,1)==1, design_bb(:,2)==1, design_bb(:,3)==1, all(design_bb==0,2)};

% Compute log power for full (1) and blank (4) epochs, at the specified frequencies
full.sl  = sl_ts(:,:,condEpochsSL{1});
blank.sl = sl_ts(:,:,condEpochsSL{4});

full.bb  = bb_ts(:,:,condEpochsBB{1});
blank.bb = bb_ts(:,:,condEpochsBB{4});

% Get SL and BB amplitudes
full.sl = log10(getstimlocked(full.sl,slFreq+1));  % Amplitude (so not squared). Square values to get units of power
blank.sl = log10(getstimlocked(blank.sl,slFreq+1)); % Amplitude (so not squared). Square values to get units of power

full.bb  = log10(getbroadband(full.bb,keepFrequencies,fs)); % Broadband data is already in units of power
blank.bb  = log10(getbroadband(blank.bb,keepFrequencies,fs)); % Broadband data is already in units of power

data = {full, blank};
assert(isequal(badChannels,badChannels0))

end