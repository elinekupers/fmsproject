function [data, badChannels] = loadData(dataDir, whichSession, type)
% Load data for analysis and figures of the manuscript:
%   A visual encoding model links magnetoencephalography
%   signals to neural synchrony in human cortex.
%       by Kupers, Benson, Winawer (YEAR) JOURNAL.
%
% INPUTS:
%   [dataDir]         : (str)  directory to get observed MEG data
%   [whichSession]    : (int)  data session number for subject 
%   [type]            : (str)  load different types of data, choice from:
%                               - 'SNR': signal-to-noise ratio
%                               - 'amplitudes': 12 Hz stimulus-locked
%                                   amplitudes & broadband power (60-150
%                                   Hz, excl harmonics)
%                               - 'amplitudesHigherHarmonics': 
%                                   higher stimulus-locked harmonics in the
%                                   same frequency range as broadband power
%                                   (72,84,96,108,132,144 Hz) & broadband
%                                   power (60-150 Hz, excl harmonics)
%                               - 'timeseries': sensor time series, so not
%                                   summarized into a stimulus-locked or
%                                   broadband data component
% OUTPUTS:
%   [data]           : (struct) requested data
%                               - 'SNR': cell with sl and bb snr
%                                       {(1xsensors),(1xsensors)}
%                               - 'amplitudes': struct with sl and bb
%                                   contrast (full-blank) (1xsensors)
%                               - 'timeseries': struct with  
%                                   ts (sensors x epochs x time),
%                                   condEpochsFull  (1xepochs)
%                                   condEpochsBlank (1xepochs)
%                                   badChannels     (1xsensors)
%   [badChannels]    : (vector) boolean marking all sensors labeled as bad
%
%
% By Eline Kupers, NYU (2017)


switch type
    case 'SNR'
        bb = load(sprintf(fullfile(dataDir, 's%02d_denoisedData_bb.mat'),whichSession));
        sl = load(sprintf(fullfile(dataDir, 's%02d_denoisedData_sl.mat'),whichSession));

        % Pick full field condition (first one, or the only one)
        if whichSession >8; whichCondition = 1; else whichCondition = [1 0 0]; end
    
        % get stimulus-locked snr
        snr_sl.coh = getsignalnoise(sl.results.finalmodel(1), whichCondition, 'SNR',sl.badChannels);
        snr_sl.coh = to157chan(snr_sl.coh,~sl.badChannels,'nans');
       
        % get broadband snr for before
        snr_bb = getsignalnoise(bb.results.finalmodel(1), whichCondition, 'SNR',bb.badChannels);
        snr_bb = to157chan(snr_bb,~bb.badChannels,'nans');
     
        data = {snr_sl.coh, snr_bb};
        
    case {'amplitudes', 'amplitudesHigherHarmonics'}
        
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
        if ~isequal(badChannels,badChannels0)  
            badChanIdx = union(find(badChannels),find(badChannels0));
            badChannels = false(size(badChannels0));
            badChannels(badChanIdx) = true;
        end

        % For some reason the preprocessing steps didn't take out bad
        % channel 98 in several sessions, so we do it manually here.
        if intersect(whichSession, [1:7,10,12]) 
            badChannels(98) = true; 
        end
        
        fprintf('(%s): Data session %d - Selected bad channels are: %s\n', mfilename, whichSession, sprintf('%d ', find(badChannels)));              
        
        % Remove bad channels and bad epochs from data and conditions
        sensorData = sensorData(:,~badEpochs, ~badChannels);

        % Permute sensorData for denoising
        sensorData = permute(sensorData, [3 1 2]);

        sl_ts = sensorData;
        bb_ts = denoisedts_bb{1};
        
        if intersect(whichSession, [1:7,10,12]) 
            bb_ts(98,:,:) = []; 
        end
        
        % there is a difference between the datasets in terms of scaling units
        % session 1-6 are in fempto Tesla  whereas 7-12 are in Tesla
        if any(intersect(whichSession, 9:14))
            assert(max(sl_ts(:), [], 'omitnan') < 1^-12)
            
            sl_ts = sl_ts .* 10^15;
            bb_ts = bb_ts .* 10^15;
        end

        design = zeros(size(conditions,1),3);
        design(conditions == 1,1) = 1; % Full
        design(conditions == 5,2) = 1; % Right
        design(conditions == 7,3) = 1; % Left

        % Define full field condition (first column (2nd and 2rd are right and
        % left, zeros are blanks)
        condEpochsFull = design(:,1);

        % Define blank epochs following fullfield epochs and remove bad epochs
        condEpochsBlank = logical([zeros(6,1); design(1:end-6,1)]);

        % Get rid of bad epochs
        condEpochsBlank = condEpochsBlank(~badEpochs,:);
        condEpochsFull  = condEpochsFull(~badEpochs,:);

        if sum(condEpochsFull==1)~=sum(condEpochsBlank==1)
            idxFull  = find(condEpochsFull==1);
            idxBlank = find(condEpochsBlank==1);
            if numel(idxBlank) > numel(idxFull)
                condEpochsBlank_truncated = idxBlank(randi([1, length(idxBlank)],1, length(idxFull)));
                condEpochsFull_truncated = idxFull;
            elseif numel(idxBlank) < numel(idxFull)
                condEpochsFull_truncated  = idxFull(randi([1, length(idxFull)],1, length(idxBlank)));
                condEpochsBlank_truncated = idxBlank;
            end
        else
            condEpochsFull_truncated = find(condEpochsFull==1);
            condEpochsBlank_truncated = find(condEpochsBlank==1);
        end

        % Select timeseries in epochs of interest
        sl.full  = sl_ts(:,:,condEpochsFull_truncated);
        sl.blank = sl_ts(:,:,condEpochsBlank_truncated);

        bb.full  = bb_ts(:,:,condEpochsFull_truncated);
        bb.blank = bb_ts(:,:,condEpochsBlank_truncated);
        
        % Compute power (BB) or amplitudes (SL) for full and blank epochs at specified frequencies
        if strcmp(type, 'amplitudesHigherHarmonics')
            harmonics12Hz = [12:12:48,72:12:150]; %72,84,96,108,132,144
            harmonicsFreqIdx = @(x) x(harmonics12Hz+1);

            % Stimulus locked using incoherent spectrum
            sl.full = sqrt(getbroadband(sl.full,harmonicsFreqIdx, fs));  % Units of power, square root to get amplitudes
            sl.full = to157chan(sl.full, ~badChannels,'nans');

            sl.blank = sqrt(getbroadband(sl.blank,harmonicsFreqIdx, fs)); % Units of power, square root to get amplitudes
            sl.blank = to157chan(sl.blank, ~badChannels,'nans');
            
        else
            % Stimulus locked using incoherent spectrum
            sl.full = getstimlocked(sl.full,slFreq+1);  % Amplitude (so not squared). Square values to get units of power
            sl.full = to157chan(sl.full, ~badChannels,'nans');

            sl.blank = getstimlocked(sl.blank,slFreq+1); % Amplitude (so not squared). Square values to get units of power
            sl.blank = to157chan(sl.blank, ~badChannels,'nans');

            % Stimulus locked using coherent spectrum
            sl.full_coherent  = getstimlocked_coherent(sl_ts,slFreq+1, condEpochsFull_truncated);  % Amplitude (so not squared). Square values to get units of power
            sl.full_coherent  = to157chan(sl.full_coherent, ~badChannels,'nans');

            sl.blank_coherent = getstimlocked_coherent(sl_ts,slFreq+1, condEpochsBlank_truncated); % Amplitude (so not squared). Square values to get units of power
            sl.blank_coherent = to157chan(sl.blank_coherent, ~badChannels,'nans');
        end
        
        % Broadband
        bb.full   = getbroadband(bb.full ,keepFrequencies,fs); % Broadband data is already in units of power
        bb.full   = to157chan(bb.full, ~badChannels,'nans');

        bb.blank  = getbroadband(bb.blank,keepFrequencies,fs); % Broadband data is already in units of power       
        bb.blank   = to157chan(bb.blank, ~badChannels,'nans');
  
        % Put in data struct
        data.sl = sl;
        data.bb = bb;
        
     case 'timeseries'
         
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
         
         % Load data
         load(sprintf(fullfile(dataDir, 's%02d_conditions.mat'),whichSession));
         load(sprintf(fullfile(dataDir, 's%02d_sensorData.mat'),whichSession));

         % preprocessing parameters (see nppPreprocessData)
         varThreshold        = [0.05 20];
         badChannelThreshold = 0.2;
         badEpochThreshold   = 0.2;
         dataChannels        = 1:157;
         
         % Preprocess raw sensordata
         [sensorData, badChannels, badEpochs] = nppPreprocessData(sensorData(:,:,dataChannels), ...
             varThreshold, badChannelThreshold, badEpochThreshold, false);
         
         % ---- Define first epochs in order to remove later ------------------
         badEpochs(1:6:end) = 1;
         
         % For some reason the preprocessing steps didn't take out bad
         % channel 98 in several sessions, so we do it manually here.
         if intersect(whichSession, [1:7,10,12]) 
            badChannels(98) = true; 
         end
         
         % Remove bad channels and bad epochs from data and conditions
         sensorData = sensorData(:,~badEpochs, ~badChannels);
         
         % Permute sensorData for denoising
         sensorData = permute(sensorData, [3 1 2]);
         
        % Define design matrix
        design = false(size(conditions,1),3);
        design(conditions == 1,1) = true; % Full
        design(conditions == 5,2) = true; % Right
        design(conditions == 7,3) = true; % Left

        % Define full field condition (first column (2nd and 2rd are right and
        % left, zeros are blanks)
        condEpochsFull = design(:,1);

        % Define blank epochs following fullfield epochs and remove bad epochs
        condEpochsBlank = logical([zeros(6,1); design(1:end-6,1)]);

        % Get rid of bad epochs
        condEpochsBlank = condEpochsBlank(~badEpochs,:);
        condEpochsFull  = condEpochsFull(~badEpochs,:);
        
        if any(intersect(whichSession, 9:14))
            assert(max(sensorData(:), [], 'omitnan') < 1^-12)
            
            sensorData = sensorData .* 10^15;
            sensorData = sensorData .* 10^15;
        end
        
        % Add timeseries and additional info to struct
        data.ts              = sensorData;
        data.condEpochsFull  = condEpochsFull;
        data.condEpochsBlank = condEpochsBlank;
        data.badChannels     = badChannels;

end