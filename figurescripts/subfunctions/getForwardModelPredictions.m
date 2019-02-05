function w = getForwardModelPredictions(G, template, f, n, nrEpochs)

% Function to compute forward model predictions from gain matrix and templates. 
% Templates will contain coherent and incoherent simulated signals. 

% INPUTS:
% G                  : [matrix] Gain matrix [nSensors x nSources]
% template           : [matrix] template to get predictions for [1xnSources]
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


%% Create cycles of a sine wave

% Create time vector
t  = (1:n)/n;   % s

% Preallocate space
signalCoherent   = zeros([size(template,2), n, nrEpochs]);
signalIncoherent = zeros(size(signalCoherent));
signalMix        = zeros(size(signalCoherent));

% Get nr of vertices in template that are actually used
nrVertices = sum(template);

% Define von Mises params
theta = 0; % mean preferred phase
kappa.incoh = 0;
kappa.coh = 20*pi;
kappa.mix = pi;

% Sample phases from three different von Mises 
phase.coh   = circ_vmrnd(theta, kappa.coh, nrEpochs);
phase.incoh = circ_vmrnd(theta, kappa.incoh, [nrEpochs, nrVertices]);
phase.mix   = circ_vmrnd(theta, kappa.mix, [nrEpochs, nrVertices]);

% Create time series, add different phases sampled from von Mises later
ts = repmat((2*pi*f*t), nrEpochs, nrVertices);
ts = reshape(ts, [nrEpochs, n, nrVertices]);

% Coherent signal contains the same phase per vertex, but differs
% slightly per epoch (depending on the width of the von Mises
% distribution)
tsCoherent = sin(ts+phase.coh);
tsCoherent = permute(tsCoherent, [3, 2, 1]);

signalCoherent(find(template),:,:) =  tsCoherent;

% Incoherent signal gets a random phase every vertex and every epoch
phase.incoh = reshape(phase.incoh, [nrEpochs, nrVertices]);

ts = permute(ts, [1,3,2]);
tsIncoherent = sin(ts + phase.incoh); 
tsIncoherent = permute(tsIncoherent, [2, 3, 1]); 

signalIncoherent(find(template),:,:) = tsIncoherent;
        
% Signal mix 
tsMix = sin(ts + phase.mix);
tsMix = permute(tsMix, [2, 3, 1]); 
signalMix(find(template),:,:) = tsMix;

    
%% Predicted sensor time series
for epoch = 1:nrEpochs
    w.c(:,:,epoch) = G*signalCoherent(:,:,epoch);
    w.i(:,:,epoch) = G*signalIncoherent(:,:,epoch);
    w.m(:,:,epoch) = G*signalMix(:,:,epoch);
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

