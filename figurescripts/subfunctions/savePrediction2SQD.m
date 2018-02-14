function ok = savePrediction2SQD(fmsRootPath, subject, prediction, fName)

% Create folder to save sqd
if ~exist(fullfile(fmsRootPath, 'data', subject),'dir')
    mkdir(fullfile(fmsRootPath, 'data', subject));
end

% Get example meg sqd file 
exampleSQD = (fullfile(fmsRootPath,'external', 'exampleFile.sqd'));

% Add additional NaNs to pad the sensor space to 192
dataToSave = reshape(prediction, size(prediction,1), []);
dataToSave = [dataToSave; NaN(192-size(dataToSave,1),size(dataToSave,2))];

% Create new sqd filename
newFile = fullfile(fmsRootPath, 'data', subject, [fName '.sqd']);

% Write padded prediction to new sqd file
sqdwrite(exampleSQD, newFile, dataToSave');

if exist(newFile, 'file') 
    ok = 1;
end

return