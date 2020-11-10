function tsmatrix = ecogTSeriesVector2TSeriesMatrix(ts, onsets)
% Convert a time series vector in a matrix of epochs x time points
%
% tsmatrix = ecogTSeriesVector2TSeriesMatrix(ts, onsets)
%
% Inputs
%  ts:          time series vector for a single run
%  onsets:      vector of sample numbers corresponding to start of epochs
%
% Outputs
%   tsmatrix :  matrix of epochs x time points


nsamples = median(diff(onsets)); % epoch length (number of temporal samples)
tsmatrix = zeros(nsamples,length(onsets));
for ii = 1:length(onsets)
    tsmatrix(:,ii) = ts((0:nsamples-1) + onsets(ii));
end

return