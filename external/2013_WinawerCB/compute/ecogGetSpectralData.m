function [Amp, Phase, f] = ecogGetSpectralData(tsmatrix, T, fmax, useHann)
% Calculate the spectrum via short-time fourier transform of an ECoG
% experiment.
%
%   [Amp, Phase, f] = ecogGetSpectralData(tsmatrix, T, fmax, useHann)
%
%
% Inputs
%  tsmatrix:    time series vector for a single run
%  fmax:        maximum temporal frequency (scalar)
%  useHann:     boolean. if true, then use a hanning window for spectra
%
% Return
%   Amp: amplitude or power spectrum (trials x frequencies)
%   f: row vector of frequencies
%   Phase: phase spectrum (trials x frequencies)
%% Set up
if notDefined('fmax'),        fmax        = 150;   end
if notDefined('useHann'),     useHann     = true;  end
       

%% Hanning window

% time series should be time points by epochs 
nt = size(tsmatrix,1); 

if useHann, H = hann(nt);
else        H = ones(nt,1); end

%% Spectrum

% frequencies of interest 
f = 0:fmax;

% The frequencies need to be adjusted to account for the display refresh
% rate. Because the rate is close to, but not precisely 60 Hz, the time
% that a stimulus aperture is displayed is close to, but not precisely 1 s.
% the actual time window (T) is usually slightly more than 1. we adjust the
% frequencies down by the same amount. hence the steady state frequency is
% not 15 Hz, but rather 15/T Hz.
f = f/T;

% put the ts matrix into a new, windowed matrix. we keep the original ts
% matrix untouched so that we can return it, if requested.
tsWindowed = bsxfun(@times, tsmatrix, H);

FT    = fft(tsWindowed)';
FT    = FT(:, 1:length(f));
Amp   = abs(FT)/nt*2;
Phase = angle(FT);

return
