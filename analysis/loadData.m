function data = loadData(dataDir, whichSession, varargin)
% Load data for analysis and figures of the manuscript:
%   A visual encoding model links magnetoencephalography
%   signals to neural synchrony in human cortex.
%       by Kupers, Benson, Winawer (YEAR) JOURNAL.
%
% INPUTS:
%   dataDir           : (str)  directory to get observed MEG data
%   whichSession      : (int)  data session number for subject
%   [type]            : (str)  load different types of data, choice from:
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
%                               - 'amplitudesCoherentSpectrum': amplitudes
%                                   using the coherent spectrum. i.e. first
%                                   taking the average across epochs in
%                                   time, then compute abs fft.
%                               - 'amplitudesBB10hz': compute broadband
%                                   response in 10 Hz bins from 60-150
%                                   excluding sl harmonics.
%  [useSLPower]      :  (bool) convert SL amplitudes to power units
%  [useBBPower]      :  (bool) convert BB amplitudes to power units
%
% OUTPUTS:
%   data             : (struct) requested data
%                               - 'amplitudes' or 'amplitudesHigherHarmonics':
%                                   struct with sl and bb amps, bootstraps
%                                   and SNR from contrast (full-blank)
%                               - 'timeseries': struct with
%                                   ts (sensors x epochs x time),
%                                   condEpochsFull  (1xepochs)
%                                   condEpochsBlank (1xepochs)
%                                   badChannels     (1xsensors)
%   badChannels      : (vector) boolean marking all sensors labeled as bad
%
%
% By Eline Kupers, NYU (2017)

p = inputParser;
p.KeepUnmatched = true;
p.addRequired('dataDir', @ischar);
p.addRequired('whichSession', @isnumeric);
p.addParameter('type','amplitudes', ...
    @(x) any(validatestring(x,{'amplitudes', 'amplitudesHigherHarmonics', ...
    'timeseries', 'amplitudesCoherentSpectrum','amplitudesBB10hz'})));
p.addParameter('useSLPower', false, @islogical)
p.addParameter('useBBPower', false, @islogical)
p.parse(dataDir, whichSession, varargin{:});

% Rename variables
type            = p.Results.type;
useSLPower      = p.Results.useSLPower;
useBBPower      = p.Results.useBBPower;

switch type
    
    case {'amplitudes', 'amplitudesHigherHarmonics','amplitudesCoherentSpectrum', 'amplitudesBB10hz'}
        
        % Number of bootstraps
        nBoot = 1000;
        
        % Define parameters to get SL and BB frequencies
        fs           = 1000;         % Sample rate
        f            = 0:150;        % Limit frequencies to [0 150] Hz
        slFreq       = 12;           % Stimulus-locked frequency
        tol          = 1.5;          % Exclude frequencies within +/- tol of sl_freq
        slDrop       = f(mod(f, slFreq) <= tol | mod(f, slFreq) > slFreq - tol);
        
        % Exclude all frequencies below 60 Hz when computing broadband power
        lfDrop       = f(f<60);
        
        if strcmp(type, 'amplitudesBB10hz')
            start_range = 60:10:140;
            for ii = 1:length(start_range)
                lfDrop = f(f<start_range(ii));
                hfDrop = f(f>start_range(ii)+10);
                
                % Define the frequenies and indices into the frequencies used to compute
                % broadband power
                [~, abIndex(ii).freq] = setdiff(f, [slDrop lfDrop hfDrop]);

            end
        else
            
            % Define the frequenies and indices into the frequencies used to compute
            % broadband power
            [~, abIndex] = setdiff(f, [slDrop lfDrop]);
            
            % Create function handles for the frequencies that we use
            bbFreqIdx    = @(x) x(abIndex);
        end
        
        %% Load data
        load(sprintf(fullfile(dataDir, 's%02d_conditions.mat'),whichSession));
        load(sprintf(fullfile(dataDir, 's%02d_sensorData.mat'),whichSession)); % not denoised
        load(sprintf(fullfile(dataDir, 's%02d_denoisedData_bb.mat'),whichSession));
        load(sprintf(fullfile(dataDir, 's%02d_denoisedts.mat'),whichSession)); % denoised
        
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
        % This channel is permanently broken and should not be used.
        if intersect(whichSession, [1:7,10,12])
            badChannels(98) = true;
        end
        
        fprintf('(%s): Data session %d - Selected bad channels are: %s\n', mfilename, whichSession, sprintf('%d ', find(badChannels)));              
        fprintf('(%s): Data session %d - Total bad channels: %d (%1.2f %%)\n', mfilename, whichSession, sum(badChannels),100*(sum(badChannels)/length(dataChannels)));              
        fprintf('(%s): Data session %d - Total bad epochs: %d (%1.2f %%)\n', mfilename, whichSession, sum(badEpochs),100*(sum(badChannels)/size(sensorData,2)));              
        
        % Remove bad channels and bad epochs from data and conditions
        sensorData = sensorData(:,~badEpochs, ~badChannels);
        
        % Permute sensorData to match dimensions of denoised data
        sensorData = permute(sensorData, [3 1 2]);
        
        sl_ts = sensorData;
        bb_ts = denoisedts_bb{1};
        
        % Remove bad channel 98 also for denoised broadband data
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
        
        % Get experimental design with stimulus used in each epoch
        design = zeros(size(conditions,1),3);
        design(conditions == 1,1) = 1; % Full
        design(conditions == 5,2) = 1; % Right
        design(conditions == 7,3) = 1; % Left
        
        % Define full field condition (first column (2nd and 3rd are right and
        % left, zeros are blanks)
        condEpochsFull = logical(design(:,1));
        
        % Define blank epochs following fullfield epochs
        condEpochsBlank = logical([zeros(6,1); design(1:end-6,1)]);
        
        % Remove bad epochs
        condEpochsFull  = condEpochsFull(~badEpochs,:);
        condEpochsBlank = condEpochsBlank(~badEpochs,:);
        
        % Select timeseries in epochs of interest
        sl.full  = sl_ts(:,:,condEpochsFull);
        sl.blank = sl_ts(:,:,condEpochsBlank);
        
        bb.full  = bb_ts(:,:,condEpochsFull);
        bb.blank = bb_ts(:,:,condEpochsBlank);
        
        %% Compute power (BB) or amplitudes (SL) for full and blank epochs at specified frequencies
        if strcmp(type, 'amplitudesHigherHarmonics')
            harmonics12Hz = [slFreq:slFreq:(4*slFreq),(6*slFreq):slFreq:(9*slFreq), (11*slFreq), (12*slFreq)]; %12,24,36,48,72,84,96,108,132,144
            slFreqIdx = @(x) x(harmonics12Hz+1);
            
        else
            slFreqIdx = slFreq+1; % 12 Hz
        end
        
        % Bootstrapping to get std of mean sl for blank and full epochs
        if strcmp(type, 'amplitudesHigherHarmonics')
            % Stimulus locked using incoherent spectrum across all 12 Hz
            % harmonics (except line noise at 60 and 120 Hz)
            % Note: we use the getbroadband function to compute the geomean
            % across multiple frequencies
            sl.amps_full  = to157chan(getbroadband(sl.full, slFreqIdx,fs),~badChannels,'nans');
            sl.amps_blank = to157chan(getbroadband(sl.blank, slFreqIdx,fs),~badChannels,'nans');  
            
        elseif strcmp(type, 'amplitudesCoherentSpectrum')
            % Stimulus locked using coherent spectrum
            % Amplitude (so not squared). Square values to get units of power
            sl.amps_full  = to157chan(getstimlocked_coherent(sl_ts, slFreqIdx, condEpochsFull),~badChannels,'nans');
            sl.amps_blank = to157chan(getstimlocked_coherent(sl_ts, slFreqIdx, condEpochsBlank),~badChannels,'nans');
            epochsFull     = find(condEpochsFull);
            epochsBlank    = find(condEpochsBlank);
            epBootFullidx  = randi(length(epochsFull),nBoot,length(epochsFull));
            epBootBlankidx = randi(length(epochsBlank),nBoot,length(epochsBlank));
            
            for ii = 1:nBoot
                sl.boot_amps_full_mn(ii,:)  = to157chan(getstimlocked_coherent(sl_ts, slFreqIdx, epochsFull(epBootFullidx(nBoot,:))),~badChannels,'nans');
                sl.boot_amps_blank_mn(ii,:) = to157chan(getstimlocked_coherent(sl_ts, slFreqIdx, epochsBlank(epBootBlankidx(nBoot,:))),~badChannels,'nans');
            end
            sl.boot_amps_diff_mn  = sl.boot_amps_full_mn - sl.boot_amps_blank_mn;
            
        else
            % Stimulus locked using incoherent spectrum.
            % Amplitude (so not squared). Square values to get units of power
            sl.amps_full = to157chan(getstimlocked(sl.full, slFreqIdx),~badChannels,'nans');
            sl.amps_blank = to157chan(getstimlocked(sl.blank, slFreqIdx),~badChannels,'nans');   
        end
        
        %% Compute sample mean of broadband response
        if strcmp(type, 'amplitudesBB10hz')
            
            for jj = 1:length(abIndex)
                % Create function handles for the frequencies that we use
                thesebbIndices = abIndex(jj).freq;
                bbFreqIdx    = @(x) x(thesebbIndices);
                % Compute BB amplitude for each epoch (no bootstrapping)
                bb.amps_full(:,:,jj)  = to157chan(getbroadband(bb.full, bbFreqIdx, fs),~badChannels,'nans');
                bb.amps_blank(:,:,jj) = to157chan(getbroadband(bb.blank, bbFreqIdx, fs),~badChannels,'nans');
            end
        else
            % Compute BB amplitude for each epoch (no bootstrapping)
            bb.amps_full  = to157chan(getbroadband(bb.full, bbFreqIdx, fs),~badChannels,'nans');
            bb.amps_blank = to157chan(getbroadband(bb.blank, bbFreqIdx, fs),~badChannels,'nans');
        end
        
        % If requested, square values to get power
        if useSLPower
            sl.amps_full = sl.amps_full.^2;
            sl.amps_blank = sl.amps_blank.^2;
        end
        
        if useBBPower
            bb.amps_full = bb.amps_full.^2;
            bb.amps_blank = bb.amps_blank.^2;
        end

        
        % Define bootstrap function for averaging across epochs
        meanCondition_fun = @(epochData) mean(epochData,1, 'omitnan');
        
        if ~strcmp(type, 'amplitudesCoherentSpectrum')
            % Bootstrap diff Full and Blank epochs for SL and BB
            % (data are in units of amplitude as well)
            sl.boot_amps_full_mn  = bootstrp(nBoot,meanCondition_fun,sl.amps_full);
            sl.boot_amps_blank_mn = bootstrp(nBoot,meanCondition_fun,sl.amps_blank);
            sl.boot_amps_diff_mn  = sl.boot_amps_full_mn - sl.boot_amps_blank_mn;
        end
        
        % Function to get difference mean full and mean blank from data without bootstrapping (sample mean)
        meanDiffFullBlank_fun = @(epochDataFull,epochDataBlank) ...
               mean(epochDataFull,1, 'omitnan') - mean(epochDataBlank,1, 'omitnan');
           
        if strcmp(type, 'amplitudesBB10hz')
           sl.amps_diff_mn = meanDiffFullBlank_fun(sl.amps_full,sl.amps_blank);
           
           % same but for broadband 10 Hz bands 
           for jj = 1:size(bb.amps_full,3)
                % Bootstrap broadband data
                bb.boot_amps_full_mn(:,:,jj)  = bootstrp(nBoot,meanCondition_fun,bb.amps_full(:,:,jj));
                bb.boot_amps_blank_mn(:,:,jj) = bootstrp(nBoot,meanCondition_fun,bb.amps_blank(:,:,jj));
                bb.boot_amps_diff_mn(:,:,jj)  = bb.boot_amps_full_mn(:,:,jj) - bb.boot_amps_blank_mn(:,:,jj);

                bb.amps_diff_mn(:,:,jj) = meanDiffFullBlank_fun(bb.amps_full(:,:,jj),bb.amps_blank(:,:,jj));
            end
        else
            % Bootstrap broadband data
            bb.boot_amps_full_mn  = bootstrp(nBoot,meanCondition_fun,bb.amps_full);
            bb.boot_amps_blank_mn = bootstrp(nBoot,meanCondition_fun,bb.amps_blank);
            bb.boot_amps_diff_mn  = bb.boot_amps_full_mn - bb.boot_amps_blank_mn;

            % Get sample mean
            sl.amps_diff_mn = meanDiffFullBlank_fun(sl.amps_full,sl.amps_blank);
            bb.amps_diff_mn = meanDiffFullBlank_fun(bb.amps_full,bb.amps_blank);
        end
        
        % Define signal as amps_diff_mn
        sl.signal = sl.amps_diff_mn;
        bb.signal = bb.amps_diff_mn;
        
        % Get noise from bootstraps (standard error)
        sl.noise = std(sl.boot_amps_diff_mn, 1);
        bb.noise = std(bb.boot_amps_diff_mn, 1);
        
        % Get SNR
        sl.snr = sl.signal ./ sl.noise;
        bb.snr = bb.signal ./ bb.noise;
        
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
        bbFreqIdx    = @(x) x(abIndex);
        
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
        fprintf('(%s): Data session %d - Bad epochs before removing 1st epoch stim period: %d (%1.2f %%)\n', mfilename, whichSession, sum(badEpochs),100*(sum(badEpochs)/size(sensorData,2)));
         
         % ---- Define first epochs in order to remove later ------------------
         badEpochs(1:6:end) = 1;
         
         fprintf('(%s): Data session %d - Selected bad channels are: %s\n', mfilename, whichSession, sprintf('%d ', find(badChannels)));              
         fprintf('(%s): Data session %d - Total bad channels: %d (%1.2f %%)\n', mfilename, whichSession, sum(badChannels),100*(sum(badChannels)/length(dataChannels)));
         fprintf('(%s): Data session %d - Total bad epochs: %d (%1.2f %%)\n', mfilename, whichSession, sum(badEpochs),100*(sum(badEpochs)/size(sensorData,2)));
       
        
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
return