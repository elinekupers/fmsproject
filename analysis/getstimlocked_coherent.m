function sl = getstimlocked_coherent(data,indexed_frequency, epochsToUse, plotNeighboringSLFreq)

% get the response amplitude at the stimulus locked frequency taking first
% the mean (therefore it is called the coherent spectrum, instead of the
% incoherent spectrum).

% INPUT
% data   : raw time series [channels x time x epochs]
% freq   : index for stimulus locked frequency (either a number or a
%          struct)
% OUTPUT
% sl     : coherhent stimulus locked time series [time x channels]

% check input 
if notDefined('which_freq')
    indexed_frequency = 13; % corresponds to 12 Hz at 1 Hz spacing
end

if notDefined('plotNeighboringSLFreq')
    plotNeighboringSLFreq = false;
end

% Take mean across epochs
data_mn = nanmean(data(:,:,epochsToUse),3);

% Fourier transform data, take amplitude of relevant frequency
spec_coh = fft(data_mn,[],2);
sl_12       = abs(squeeze(spec_coh(:,indexed_frequency)))' / size(data,2)*2;
sl_11       = abs(squeeze(spec_coh(:,indexed_frequency-1)))' / size(data,2)*2;
sl_13       = abs(squeeze(spec_coh(:,indexed_frequency+1)))' / size(data,2)*2;

if plotNeighboringSLFreq
    sl.sl11 = sl_11;
    sl.sl12 = sl_12;
    sl.sl13 = sl_13;
else
    sl = sl_12;
end

return