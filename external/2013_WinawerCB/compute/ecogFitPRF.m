function params = ecogFitPRF(data, stimulus, useExp, xval)
% Create inputs needed for call to fitprf, and then call fitprf.
%
%   params = ecogFitPRF(data, stimulus, useExp, xval)
%
%
% Inputs:
%   data: ECoG data (epochs x 1)
%   stimulus: stimulus time series (n time points x n pixels [note that 2d image is flattened to 1D)
%   useExp: if true, use CSS pRF model otherwise fit linear pRF model
%   xval: if true, apply leave-one-run out cross validation scheme
%       (assuming 3 types of stimulus, we solve the model 3 times, and
%       report goodness of fit by concatenating the 3 left-out time series
%       and predicted time series)
%
% Outputs:
%   params structure, including a variety of fields:
%       - params, 1x6 (no xval) or nx6 (xval), with fields x,y, sigma, gain, exponent (x,y,sigma in pixels)
%       - paramsse: standard error on params (applies only to xval)
%       - r: goodness of fit measure concatenated across runs (presumably variance explained)
%       - rrun: same as r, but applied to each run separately
%       - polyparams: polynomial detreneding (typically not used for ECoG, but useful for fMRI)
%       - polymeans: something to do with polynomial detrending (not applicable for ECoG)
%       - numiters: number of iterations used by solver before arriving at solution
%       - hrf: hemodynamic response function (not applicable here; should be same as gain, i.e., 4th column in params field)
%       - betas: not applicable
%       - signal: signal predicted by the model
%       - drift: not applicable (predicted polynomial fit)

if notDefined('useExp'),            useExp = true;           end
if notDefined('xval'),              xval = false;            end


%% Define fitprf inputs

% prfmodel
x = fpDefaultInputs(stimulus);

prfmodel = x.prfmodel; % 2-d gaussian (with static non-linearity)

% hrfmodel
hrfmodel = 1; % impulse response function

% flag
flag = 1; % NA

% maxpolydeg
blanks = sum(catcell(1, stimulus),2) < eps;
if sum(blanks) > 0, maxpolydeg = NaN;    % no detrending
else                maxpolydeg = 0;  end % model the mean

% ar
ar = []; % NA

% mode
if useExp,  mode = {1 1 [0; Inf]};
else        mode = 0;  end

% maxiter
maxiter = 2000;

% wantresample
wantresample = 0; % for cross-validation

% tol
tol = x.tol;

% extraopt
extraopt = x.extraopt;
extraopt{4} = 'final';%'iter';
% extraregressors
extraregressors = [];

% hrfnormfun
hrfnormfun = [];

% derivemode
derivemode = [];

% metric
if sum(blanks) > 0, metric = @(x,y)calccod(x,y,[],[],0); % if we know the baseline
else                metric = @(x,y)calccod(x,y,[],[],1); end

% outputfcn
outputfcn = [];

%% SOLVE!!


% if cross-validate, then we do n-fold leave-one out
if xval, wantresample = 'n-fold'; end

% run the prf fit!

params = fitprf(stimulus,data,prfmodel,hrfmodel,flag,maxpolydeg,ar,mode,maxiter, ...
    wantresample,tol,extraopt,extraregressors,hrfnormfun,derivemode,metric,outputfcn);



