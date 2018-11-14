function [data, badChannels] = loadData(dataDir, whichSession, type)


switch type
    case 'SNR'
        bb = load(sprintf(fullfile(dataDir, 's%02d_denoisedData_bb.mat'),whichSession));
        sl = load(sprintf(fullfile(dataDir, 's%02d_denoisedData_sl.mat'),whichSession));

        data = {bb, sl};
        
    case 'amplitudes'
        
        % Define parameters to get SL and BB data
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
        load(sprintf(fullfile(dataDir, 's%02d_conditions.mat'),whichSession));
        load(sprintf(fullfile(dataDir, 's%02d_sensorData.mat'),whichSession));
        load(sprintf(fullfile(dataDir, 's%02d_denoisedData_bb.mat'),whichSession));
        load(sprintf(fullfile(dataDir, 's%02d_denoisedts.mat'),whichSession));


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

        % Make sure bad epochs and bad sensors are the same across two datasets
        assert(isequal(badEpochs,badEpochs0))
%         assert(isequal(badChannels,badChannels0))

        % Remove bad channels and bad epochs from data and conditions
        sensorData = sensorData(:,~badEpochs, ~badChannels);

        % Permute sensorData for denoising
        sensorData = permute(sensorData, [3 1 2]);

        sl_ts = sensorData;
        bb_ts = denoisedts_bb{1};

        design = zeros(size(conditions,1),3);
        design(conditions == 1,1) = 1; % Full
        design(conditions == 5,2) = 1; % Right
        design(conditions == 7,3) = 1; % Left

        % Define full field condition (first column (2nd and 2rd are right and
        % left, zeros are blanks)
        condEpochsFull = design(:,1)==1;

        % Define blank epochs following fullfield epochs and remove bad epochs
        condEpochsBlank = logical([zeros(6,1); design(1:end-6,1)]);

        % Get rid of bad epochs
        condEpochsBlank = condEpochsBlank(~badEpochs,:);
        condEpochsFull  = condEpochsFull(~badEpochs,:);

        if sum(condEpochsFull==1)~=sum(condEpochsBlank==1)
            if sum(condEpochsBlank==1) > sum(condEpochsFull==1)
                idx = find(condEpochsBlank);
                idx = idx(1:length(find(condEpochsFull)));
                condEpochsBlank = condEpochsBlank(idx);
            else
                idx = find(condEpochsFull);
                idx = idx(1:sum(condEpochsBlank));
                condEpochsFull = condEpochsFull(idx);
            end
        end

        % Select timeseries in epochs of interest
        sl.full  = sl_ts(:,:,condEpochsFull);
        sl.blank = sl_ts(:,:,condEpochsBlank);

        bb.full  = bb_ts(:,:,condEpochsFull);
        bb.blank = bb_ts(:,:,condEpochsBlank);

        % Compute log power for full and blank epochs at specified frequencies
        sl.full = getstimlocked(sl.full,slFreq+1);  % Amplitude (so not squared). Square values to get units of power
        sl.blank = getstimlocked(sl.blank,slFreq+1); % Amplitude (so not squared). Square values to get units of power

        sl.full_coherent = getstimlocked_coherent(sl_ts,slFreq+1, condEpochsFull);  % Amplitude (so not squared). Square values to get units of power
        sl.blank_coherent = getstimlocked_coherent(sl_ts,slFreq+1, condEpochsBlank); % Amplitude (so not squared). Square values to get units of power
        
        bb.full  = getbroadband(bb.full ,keepFrequencies,fs); % Broadband data is already in units of power
        bb.blank  = getbroadband(bb.blank,keepFrequencies,fs); % Broadband data is already in units of power

        data.sl = sl;
        data.bb = bb;

end