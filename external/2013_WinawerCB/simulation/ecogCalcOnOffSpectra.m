function [on, off] = ecogCalcOnOffSpectra(on, off, useHann, calcPower)
% Calculate spectra from time series for on-off experiments
% 
%   [on, off] = ecogCalcOnOffSpectra(on, off, useHann, calcPower)
%
% Inputs
%   on: struct containing a field called 'signal', which is a matrix 
%       num time points x epochs
%   off: same, but for off epochs
%   useHann: (boolean) If true, window each epoch with a hann function
%   calcPower: (boolean) If true, return spectra as squared amplitude
%                       rather than amplitude
% Outputs
%   on:  struct containing time series and spectra for on epochs
%   off: struct containing time series and spectra for off epochs
%
% Example
%   t = .001:.001:1;
%   ntrials = 50;
%   on.signal  = randn(length(t),ntrials) + sin(repmat(t',1, ntrials)*2*pi*10);
%   off.signal = randn(length(t),ntrials) ;
%   [on off] = ecogCalcOnOffSpectra(on, off);
%

if notDefined('useHann'),   useHann = false; end
if notDefined('calcPower'), calcPower = false; end
nt = size(on.signal, 1);

if useHann,
    H  = hann(size(on.signal, 1));    
    on.signalWindowed = bsxfun(@times, on.signal, H);
    off.signalWindowed = bsxfun(@times, off.signal, H);
else
    on.signalWindowed = on.signal;
    off.signalWindowed = off.signal;
end
on.spectra  = abs(fft(on.signal, [], 1))/nt*2;
off.spectra = abs(fft(off.signal,[], 1))/nt*2;

if calcPower,
    on.spectra = on.spectra.^2;
    off.spectra = off.spectra.^2;
end

% on.meanFFT  = exp(mean(log(on.spectra), 2));
% off.meanFFT = exp(mean(log(off.spectra), 2));
on.meanFFT  = mean(on.spectra, 2);
off.meanFFT = mean(off.spectra, 2);

on.fftMean  = abs(fft(mean(on.signalWindowed,2)))/nt*2;
off.fftMean = abs(fft(mean(off.signalWindowed,2)))/nt*2;

if calcPower,
    on.fftMean = on.fftMean.^2;
    off.fftMean = off.fftMean.^2;
end


on.meanTs   = mean(on.signal, 2);
off.meanTs  = mean(off.signal, 2);

