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

signalCoherent  = zeros([size(template,2), n, nrEpochs]);
signalIncoherent = zeros(size(signalCoherent));

for epoch = 1:nrEpochs
    
    % Coherent signal gets a fixed phase for each vertex, in one epoch (but
    % different across epochs)
    thisEpochCoherentPhase = rand*2*pi;
    signalCoherent(:,:,epoch) = template' * sin(2*pi*f*t + thisEpochCoherentPhase);
    
    % Incoherent signal gets a random phase every vertex and every epoch
    for ii = find(template)
        signalIncoherent(ii,:,epoch) = sin(2*pi*f*t + rand*2*pi);
    end
    
    %% Predicted sensor time series
    w.c(:,:,epoch) = G*signalCoherent(:,:,epoch);
    w.i(:,:,epoch) = G*signalIncoherent(:,:,epoch);
    
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

