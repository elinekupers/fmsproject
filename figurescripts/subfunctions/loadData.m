function data = loadData(dataDir, whichSubject)

bb = load(sprintf(fullfile(dataDir, 's%02d_denoisedData_bb.mat'),whichSubject));
sl = load(sprintf(fullfile(dataDir, 's%02d_denoisedData_sl.mat'),whichSubject));
load(sprintf(fullfile(dataDir, 's%02d_conditions.mat'),whichSubject));
data = {bb,sl};

end