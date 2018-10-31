function preprocessMEGData(whichSubject)

% This function is to preprocess the SSMEG data:
% 1. Load raw sqd
% 2. Get triggers
% 3. Make epochs with the given triggers
% 4. Save sensorData and conditions in separate mat files

% INPUTS:
% whichSubject      : [integer] Subject dataset you want to load.
%                     The fullOnly subfolder has subjects 1-4.
%                     The general SSMEG folder (set subfolder to empty) 
%                     has subjects 1-8.

% OUTPUTS:
% None

% % Example: 
% preprocessMEGData(4)

%% Set analysis variables
projectPth                    = '/Volumes/server/Projects/MEG/SSMEG/';
subFolder                     = 'fullOnly'; 

saveTimeSeries                = true;        % Save epoched time series?
triggerChannels               = 161:164;
dataChannels                  = 1:157;

fs                            = 1000;        % sample rate (Hz)
epochTime                     = [0 1];        % start and end of epoch, relative to trigger, in seconds

%% To run script, you need the Field trip toolbox

% Add fieldtrip path
% meg_add_fieldtrip_paths('/Volumes/server/Projects/MEG/code/fieldtrip', 'yokogawa_defaults');

% Find subjects for this project
subjectPths = dir(fullfile(projectPth, subFolder, '*SSMEG_*'));

%% -------------------------------------
% ------------ Load in data ------------
% --------------------------------------

dataPth = fullfile(projectPth, subFolder, subjectPths(whichSubject).name, 'raw');

[ts, megfiles] = meg_load_sqd_data(dataPth, '*SSMEG*');
hdr = ft_read_header(fullfile(megfiles(1).folder, megfiles(1).name));

%% -------------------------------------
% ------------ Get triggers ------------
% --------------------------------------

if (~strcmp(subFolder,'fullOnly')) && (whichSubjects == 5) 
    pdChan  = 192; % Photodiode channel
    numRuns = 6;  % Per single flicker/blank period
%     num_epoch_time_pts = 1000;
    
    % This function is specifically made for session 8, look inside the
    % code if you want to use it for a different session!
    [trigger] = ssmeg_get_triggers_from_photodiode(pdChan, fs, numRuns, ts);    
else
      
    trigger = meg_fix_triggers(ts(triggerChannels,:)');
    
end

if strcmp(subFolder,'fullOnly')
    if whichSubject == 4
        
        startOfRuns = find(diff(find(trigger))>85); % 85 because of 1/12 flicker
        inds = find(trigger);

        % add lost trigger in the middle of run 9
        trigger(inds(startOfRuns(9))+82)=2; %
        
        % remove eighth run with noise
        trigger(inds(startOfRuns(7)+1):inds(startOfRuns(8))) = 0;
        
    end
    
    onsets = ssmeg_trigger_2_onsets(trigger, 7, 'meg');
else
    onsets = ssmeg_trigger_2_onsets(trigger, whichSubject, 'meg');
end

%% -------------------------------------
% ------------ Make epochs -------------
% --------------------------------------

[sensorData, conditions] = meg_make_epochs(ts(dataChannels,:)', onsets, epochTime, fs); %#ok<ASGLU>

%% -------------------------------------
% ------------- Save data --------------
% --------------------------------------

if saveTimeSeries
    save(fullfile(projectPth, subFolder, subjectPths(whichSubject).name, 'processed', sprintf('s%02d_sensorData', whichSubject)), 'sensorData');
    save(fullfile(projectPth, subFolder, subjectPths(whichSubject).name, 'processed', sprintf('s%02d_conditions', whichSubject)), 'conditions');    
end

return




