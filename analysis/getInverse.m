function s = getInverse(fName, prediction, keep_sensors)

% Load Brainstorm inverse model
inverseModel  = load(fName);

% Reshape forward model prediction
predictionTimeseries = reshape(prediction, size(prediction,1), []);

% Create source response for each timepoint
sourceTimeseries = inverseModel.ImagingKernel*predictionTimeseries(keep_sensors,:);

% Reshape back into epoched source time series
s = reshape(sourceTimeseries, [size(sourceTimeseries,1), size(prediction,2), size(prediction,3)]);
