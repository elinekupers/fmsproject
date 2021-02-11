function [amps, phase, freq] = ecogSummarizeSpectra(Amps, Ph, T, fmax)
% Convert a matrix of fourier coefficients (frequency x epoch) into 
% vectors of ECoG compoennts (stimulus-locked and broadband)
%
% [amps, phase, freq] = ecogSummarizeSpectra(Amps, Ph, T, fmax)
%
% Inputs:
%   Amps: Matrix of Fourier amplitudes, frequencies x epoch
%   Ph:   Matrix of Fourier phases, frequencies x epoch
%   T:    Epoch length, in seconds
%   fmax: maximum frequency used to compute broadband ECoG  
%
% Outputs
%   amps: structure comprised of two vectors, ab (asyncrhonous broadband)
%           and sl (stimulus-locked), each (num epochs) x 1
%   phase: same as amps, but phase rather than amplitude. note that
%           phase.ab is an empty vector since the phase has no specific
%           meaning for the broadband signal
%   freq:  a structure with fields describing which frequencies are used
%           for various computations. See ecogGetSLandABfrequencies for
%           more detail.
%% CHECK INPUTS

% by default, use frequencies up to 150 Hz for broadband calculation
if notDefined('fmax'), fmax = 150; end

% if no Ph, make a matrix of NaN to avoid errors
if notDefined('Ph'), Ph = NaN(size(Amps)); end

% Get the stimulus-locked and broadband frequencies
freq = ecogGetSLandABfrequencies((0:fmax)/T, T);

%% Asynchronous Broadband linear fit

% We assume a power law relationship between amplitude and frequency
% (amplitude is approximatley f^n, with n usually negative) as one compoent
% of the ECoG response. Because it is a power law, the relationship is
% linear in log-log space, meaning that log(amplitude) = n * log(frequency)
% + k. Because it is a power law, the linear relationship in log-log space
% holds whether the signal is defined as amplitude or squared amplitude
% (power). The only difference is that the exponent is doubled for power
% compared to amplitude. We fit the line in log-log space rather than
% fitting the power law to the untransformed spectrum because the variance
% in the log spectrum is approximately equal across frequencies.
%
% The fitting strategy is to take the power from many epochs, and fit a
% single linear function in log-log space to determine the slope (the power
% law exponent). We then find the intercept for individual epochs to get a
% measure of the broadband power from each epoch. Only a subset of
% frequencies are used for the fits; we exclude the frequencies near the
% even harmonics of the stimulus locked frequency (typicaly 15 Hz=2f) and
% near line noise and harmonics, as well as frequencies above fmax
% (typcially 150 Hz) or below fmin (typically 8 Hz).


% store the amplitude of the asynchronus broadband (AB) and stimulus locked
% (SL) time series (we use 'amplitude loosely: in fact the AB response is
% computed based on spectral power and the SL is computed based on spectral
% amplitude)
amps    = struct('ab', [], 'sl', []); 

% store the phases (relevant only for stimlus-locked responses)
phase   = struct('ab', [], 'sl', []); 

% fit the spectral data with a line in log-log space

% Take the log of the frequencies and squared amplitude to get log-log space
Y = log(Amps(:,freq.ab_i).^2); 
X = repmat(log(freq.ab), size(Y, 1), 1);

% fit ONE line to all the concatenated trials to get a single slope
p = polyfit(X(:), Y(:), 1);

% evaluate the fitted broadband line at the stimulus-locked frequency
ssa = polyval(p, log(freq.sl));

% calculate the residuals from the linear prediction to get AB reponses on
% individual epochs
residual = bsxfun(@minus, Y, polyval(p, log(freq.ab)));

% take the mean of the residuals from each epoch and add the interecept
% from the linear fit to get the  intercept for that epoch
intercepts = mean(residual,2) + ssa;

% broadband based on linear fit (we exponentiate to have units of power,
% rather than log power)
amps.ab = exp(intercepts);

% because the broadband is comprised of many frequencies, the
% phase has no specific meaning.
phase.ab = [];


%% Stimulus-locked time series (1st harmonic only)
amps.sl  = Amps(:,freq.sl_i);
phase.sl = Ph(:,freq.sl_i);

% Subtract the AB from the Stimulus-Locked. We sqrt the AB because it was
% computed as power, whereas SL is computed as amplitude
amps.sl = amps.sl - sqrt(amps.ab);

return