function [tseriesOut, blanks] = ecogSubtractBlanks(tseriesIn, runType)
% Set the time series to have a mean of 0 for epochs with stimulus blanks
% (i.e., zero contrast)
%
% [tseriesOut, blanks] = ecogSubtractBlanks(tseriesIn, runType)
%
% Input:
%   tseriesIn: tseries matrix, epochs x runs
%   runTYpe: cell array with names of experiments, equal in length to the
%               number of columns in tseriesIn
% Outputs:
%   tseriesOut: same as tSeries in, but with a zero mean for epochs with
%               stimulus blanks
%   blanks: binary matrix indicating epochs with blanks, same size as
%               tseriesIn and tseriesOut


% Get the stimulus apertures (images x epochs, where images have been
% vectorized by stacking the columns)
stimulus = ecogGetStimulus(runType);

% identify blanks as epochs in which the mean stimulus is 0
blanks = cellfun(@(x) mean(x, 2) < 10e-6, stimulus, 'UniformOutput', false);
blanks = catcell(2, blanks);

if sum(blanks(:)) == 0,
    % do nothing - there are no blanks    
else       
    % subtract the baseline
    baseline   = mean(tseriesIn(blanks(:)));
    tseriesOut = tseriesIn - baseline;
    
end

return