function spectra = ecogGetSpectra(ts, onsets, T, fmax, useHann)
% Wrapper function to compute the fourier transform for each epoch within
% an ECoG experiment.
% 
% Inputs:
%   ts:     structure containing time series for one or more experiments,
%           with each experiment corresponding to a field in the structure.
%           Each experiment in turn is a cell array, with each cell
%           corresponding to one iteration of the experiment. Each cell is
%           1 x num time points.
%
%   onsets: a structure with the same fields as ts. Each field is comprised
%           of a cell array, with each cell corresponding to one iteration
%           of the experiment
%
%   T:      Length of one epoch (seconds)
%
%   fmax:   Store temporal frequencies up to fmax
%
%   useHann:Boolean to indicate whether to use a hanning window 
%
% Outputs
%`  spectra: matrix of Fourier amplitudes, epochs x frequencies
%


% spectra = ecogGetSpectra(ts, onsets, T, fmax, useHann)

% Names of of experiments
runTypes = fields(ts);      

% Number of different experiments (there may be repeated runs of each type)
numstim  = numel(runTypes);

% We will compute the spectra for each epoch within a run by calling
% ecogGetSpectralData and write this into spectraRuns{runnum}. We will
% then concatenate epochs across repeated runs for a given experiment type
% into spectraStimuli{stimnum}. Finally we will concatenate the spectra
% across different stimuli into a single matrix, spectra, which will
% consist of epochs x frequencies.
spectraStimuli = cell(1, numstim); 

% loop across stimulus types (i.e., experiments)
for stim = 1:numstim
    thisstim = runTypes{stim};
    numruns  = size(ts.(thisstim), 2); % number of runs per stimluus type
    spectraRuns = cell(1, numruns);
    
    % loop across repeated runs of a given stimulus
    for run = 1:numruns
        tsmatrix = ecogTSeriesVector2TSeriesMatrix(ts.(thisstim){run},...
            onsets.(thisstim){run});
        spectraRuns{run} = ecogGetSpectralData(tsmatrix, T, fmax, useHann);
    end
    
    % concatenate epochs across repeated runs
    spectraStimuli{stim} = catcell(1, spectraRuns);
end

% concatenate epochs across different experiments
spectra = catcell(1, spectraStimuli);

return