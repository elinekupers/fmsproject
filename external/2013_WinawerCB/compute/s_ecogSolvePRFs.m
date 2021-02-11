% s_ecogSolvePRFs
%
% Script to re-solve population receptive field models for ECoG electrodes
% in 4 subjects. Pre-solved pRF models are stored in the data directory,
% under the the matlab files pRF_*.mat.
% 
% The main inputs to the pRF solver are
%   * a time series vector, comprising a single value for each stimulus epoch 
%       - either the stimulus-locked amplitude or the broadband amplitude for that epoch (pre-computed)
%   * a times series of binary image masks (n time points x num pixels) 
%   * boolean to specify whether to use an exponent (CSS model) or assume a linear model
%   * boolean to specify whether to cross-validate or not
%
% Before solving, we load a data file containing the time series for all
% channels in all subjects for all stimuli and all data types. For example,
%
%   %  Load the data
%   a = load(datafile)
%
%   % Show the amplitude of the asynchronous broadband for subject 2,
%   % channel 66, for the thin bar experiment
%   figure; plot(a.amps.subj2.chan66.ab.barThin);
%
%   % Note that this data is 96x1, which means 96 distinct 1-s epochs for
%   % this experiment. To see an example of how the stimulus-locked and
%   % asynchronous broadband time series were computed from the raw time
%   % series, see the script figure2cd_AB_and_SL_timeSeries.
%
%  We then compute several kinds of models within loops
%  Loop 1: cross validate (for figure 4) or not (for figure 3)
%   Loop 2: use linear or CSS model (if no x-val, then only CSS model)
%    Loop 3: stimulus-locked or broadband component of the ECoG data
%     Loop 4: subjects 1:4
%      Loop 5: multiple channels within each subject
%  Inside the last loop, we call the fit function, ecogFitPRF
%
% Results of all models get stored in savepth

%% set up
savepth   = fullfile(ecogPRFrootPath, 'data');
datafile  = fullfile(ecogPRFrootPath, 'data', 'tseriesForPRFs.mat');

load(datafile);

% This warning will appear hundres of time unless we quiet Matlab
warning off optimlib:lsqncommon:SwitchToLargeScale

% Loop 1: Cross-validate or not (save in separate files)
for xval = [false true]
    
    % if we are not cross-validating, we always use the CSS model (allow
    % for compressive nonlinearity). if we are cross-validating, then we
    % try both with and without the nonlinearity to compare model fits.
    if ~xval, useExp = true; else useExp = [false true]; end
    
    % Loop 2: Use compressive nonlinearity or not (save in separate files)
    for ue = useExp
        
        % Loop 3: Broadband or stimlus-locked ECoG compoennt (save in
        % separate files)
        datatypes = {'ab' 'sl'};
        for dt = 1:2
            
            % initialize a struct to store results for all subjects and channels
            params   = struct('subj', [], 'chan', [], 'roi', [], 'isV1V2V3', [], ...
                'x', [], 'y', [], 's', [], 'n', [], ...
                'params', [], 'paramsse', [], 'r', [], 'rrun', [], 'polyparams', [],...
                'polymeans', [], 'numiters', [], 'hrf', [], 'betas', [],...
                'signal', [], 'drift', [], 'resp', []);
            
            if strcmp(datatypes{dt}, 'ab'), subjects = 1:4; else subjects = 1:3; end
            
            % Loop 4: Subjects (results will be concatenated in a single file)
            for  su = subjects;
                
                % get the channels for this subject
                channelnames = fieldnames(amps.(sprintf('subj%d', su)));
                channels     = cellfun(@(x) sscanf(x, 'chan%d'), channelnames);                                
                labels       = rois.(sprintf('subj%d', su));
                v1v2v3       = isV1V2V3.(sprintf('subj%d', su));
                
                % get the stimulus description
                stimTypes = fieldnames(amps.(sprintf('subj%d', su)).(channelnames{1}).ab);
                stimulus  = ecogGetStimulus(stimTypes(1:3));
                
                % Loop 5: Channels (results will be concatenated in a single file)
                for ch = 1:length(channels)
                    
                    fprintf('\n\n----------------------------------------------\n')
                    fprintf('Solving models for subject %d, channel %d (%d of %d)\n', ...
                        su, channels(ch), ch, length(channels));
                    fprintf('----------------------------------------------\n')
                    
                    % Get the data for this channel               
                    data = cell(1,3);
                    for ii = 1:3
                        data{ii} = amps.(sprintf('subj%d', su)).(channelnames{ch}).(datatypes{dt}).(stimTypes{ii});
                    end
                    
                    % *****************************************************                    
                    % *****************************************************
                    % ** FIT THE PRF MODEL. THIS IS THE MAIN COMPUTATION **
                    tmp = ecogFitPRF(data, stimulus, ue, xval); % *********
                    % *****************************************************
                    % *****************************************************

                    
                    % pull some useful numbers out of the solution
                    x0 = tmp.params(1); % prf x-center in pixels
                    y0 = tmp.params(2); % prf y-center in pixels
                    n  = tmp.params(5); % prf exponent
                    s0 = tmp.params(3)/sqrt(n); % prf size
                    
                    % convert pixels to degrees of visual angle
                    [x, y, s] = ecogPix2Deg(su, x0, y0, s0);
                    
                    % concatenate the struct
                    params.subj         = cat(1, params.subj, su);
                    params.chan         = cat(1, params.chan, channels(ch));
                    params.roi          = [params.roi labels(ch)];
                    params.isV1V2V3     = cat(1, params.isV1V2V3, v1v2v3(ch));
                    params.x            = cat(1, params.x, x);
                    params.y            = cat(1, params.y, y);
                    params.s            = cat(1, params.s, s);
                    params.n            = cat(1, params.n, n);
                    params.params       = cat(1, params.params, tmp.params);
                    params.paramsse     = cat(1, params.paramsse, tmp.paramsse);
                    params.r            = cat(1, params.r, tmp.r);
                    params.rrun         = cat(1, params.rrun, tmp.rrun);
                    params.polyparams   = cat(1, params.polyparams, tmp.polyparams);
                    params.polymeans    = cat(1, params.polymeans, tmp.polymeans);
                    params.numiters     = cat(1, params.numiters, tmp.numiters);
                    params.hrf          = cat(1, params.hrf, tmp.hrf);
                    params.betas        = cat(1, params.betas, tmp.betas);
                    params.signal       = catmatrix(2, params.signal, tmp.signal);
                    params.drift        = catmatrix(2, params.drift, tmp.drift);
                    params.resp         = catmatrix(2, params.resp, catcell(1,data));
                    
                    
                    % (re-)save with each new channel, in case of crash
                    if ue, model = 'exp'; else model = 'noexp'; end
                    if xval, model = ['xval_' model]; end
                    fname = sprintf('PRF_%s_%s', model, datatypes{dt});
                    save(fullfile(savepth, fname), 'params');
                end
                 
            end
        end
        
    end
end

% Clean up by turning the warnings back on
warning on optimlib:lsqncommon:SwitchToLargeScale