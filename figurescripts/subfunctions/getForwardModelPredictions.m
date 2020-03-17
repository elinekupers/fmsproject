function w = getForwardModelPredictions(G, template, f, n, nrEpochs, theta, kappa)

% Function to compute forward model predictions from gain matrix and templates. 
% Templates will contain coherent and incoherent simulated signals. 

% INPUTS:
% G                  : [matrix] Gain matrix [nSensors x nSources]
% template           : [matrix] Anatomical template of V1-V3 retinotopy [1xnSources]
% f                  : [int] frequency to use for simulate sine wave (Hz)
% n                  : [int] number of time points for simulated signal (ms)
% nrEpochs           : [int] number of epochs containing a simulated signal
%                       with n timepoints

% OUTPUTS:
% w                 :  [struct] forward model prediction for each template,
%                        using a coherent or incoherent signal

%% Check inputs
if ~exist('f', 'var') || isempty(f)
    f = 1;
end

if ~exist('n', 'var') || isempty(n)
    n = 1;
end

if ~exist('nrEpochs', 'var') || isempty(nrEpochs)
    nrEpochs = 1;
end

if ~exist('theta', 'var') || isempty(theta)
    theta = 0;  % mean preferred phase of von Mises (same across all three)
end

if ~exist('kappa', 'var') || isempty(kappa)
    kappa.asyn = 0; % (inverse) width of distribution
    kappa.syn = 100*pi;
    kappa.mix = pi;
end

%% Create cycles of a sine wave

% Create time vector
t  = (1:n)/n;   % s

% Preallocate space
signalSynchronous   = zeros([size(template,2), n, nrEpochs]);
signalAsynchronous  = zeros(size(signalSynchronous));
signalMix           = zeros(size(signalSynchronous));

% Get nr of vertices in template that are actually used
nrVertices = sum(template);

% Sample phases from three different von Mises 
phase.syn   = circ_vmrnd(theta, kappa.syn, nrEpochs);
phase.asyn  = circ_vmrnd(theta, kappa.asyn, [nrEpochs, nrVertices]);
phase.mix   = circ_vmrnd(theta, kappa.mix, [nrEpochs, nrVertices]);

% Create time series, add different phases sampled from von Mises later
timepoints = 2*pi*f*t;
ts = repmat(timepoints', [1 nrEpochs]);

%% Synchronous signals 
% Signals contain the same phase across vertices.
% Depending on the width of the von Mises distribution, they might slightly
% differ between per epoch

% Create synchronous sine waves
tsSynchronous = sin(ts+phase.syn');

% Add singleton dimension
tsSynchronous = repmat(tsSynchronous, [1 1 nrVertices]);

% Change order of array dimensions to correspond to preallocated array
tsSynchronous = permute(tsSynchronous, [3, 1, 2]);

% Add timeseries to selected vertices from anatomical template
signalSynchronous(find(template),:,:) =  tsSynchronous; %#ok<FNDSB>

%% Asynchronous signals
% Signals contain sine waves with random (uniform distributed) phases, for
% every vertex and every epoch

% Create more time points first
ts = repmat(timepoints', [1 nrEpochs, nrVertices]);

% Add singleton dimension
phase.asyn = reshape(phase.asyn, [1, nrEpochs, nrVertices]);

% Create asynchronous sine waves
tsAsynchronous = sin(ts+phase.asyn);

% Change order of array dimensions to correspond to preallocated array
tsAsynchronous = permute(tsAsynchronous, [3, 1, 2]); 

% Add timeseries to selected vertices from anatomical template
signalAsynchronous(find(template),:,:) = tsAsynchronous; %#ok<FNDSB>
        
%% Signal mix
% Same as for asynchronous, but now signals contains sine waves with 
% a mixture of phases biased towards a particular value, which can differ 
% for every vertex and every epoch

% Add singleton dimension
phase.mix = reshape(phase.mix, [1, nrEpochs, nrVertices]);

% Create mixture of more/less synchronous sine waves
tsMix = sin(ts+phase.mix);

% Change order of array dimensions to correspond to preallocated array
tsMix = permute(tsMix, [3, 1, 2]); 

% Add timeseries to selected vertices from anatomical template
signalMix(find(template),:,:) = tsMix; %#ok<FNDSB>

    
%% Predicted sensor time series by multiplying source activity with Gain

for epoch = 1:nrEpochs
    w.c(:,:,epoch) = G*signalSynchronous(:,:,epoch);  % synchronous or coherent
    w.i(:,:,epoch) = G*signalAsynchronous(:,:,epoch); % asynchronous or incoherent
    w.m(:,:,epoch) = G*signalMix(:,:,epoch);          % mix
end

    
return

%% Note: other way simulating coherent/incoherent signal (w/ complex numbers)

% Compute the sensor weights, w, from V1-V3 using the contrained gain field (forward model)
%     w = G_constrained*V1template'; %  Nsensors x 1;
    
%     phAmp2complex = @(r,th) r .* exp(1i*th);
  
%     template.V1StimEccenPhaseCoherent     = template.V1StimEccenAmplitudes * 0;    
%     template.V1StimEccenPhaseIncoherent   = template.V1StimEccenAmplitudes .* (rand(iter,size(template.V1StimEccenAmplitudes,2)) * 2*pi);    
    
    % Make a complex number with amplitudes and phases:
%     template.V1coherent = phAmp2complex(template.V1StimEccenAmplitudes,template.V1StimEccenPhaseCoherent);     
%     template.V1incoherent = phAmp2complex(repmat(template.V1StimEccenAmplitudes,[iter,1]),template.V1StimEccenPhaseIncoherent);

