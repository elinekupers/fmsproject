%% s_saveSummarySQD.m

% script to save several data/model outputs as an .sqd file

% Define subjects
subject         = {'wlsubj002', ... % S1 - Full, Left, Right stim experiment
                   'wlsubj004', ... % S2 - Full, Left, Right stim experiment
                   'wlsubj005', ... % S3 - Full, Left, Right stim experiment
                   'wlsubj006', ... % S4 - Full, Left, Right stim experiment
                   'wlsubj010', ... % S5 - Full, Left, Right stim experiment
                   'wlsubj011', ... % S6 - Full, Left, Right stim experiment
                   'wlsubj048', ... % S7 - Full  stim only experiment
                   'wlsubj046', ... % S8 - Full  stim only experiment
                   'wlsubj039', ... % S9 - Full  stim only experiment
                   'wlsubj059', ... % S10 - Full  stim only experiment
                   'wlsubj067', ... % S11 - Full  stim only experiment
                   'wlsubj070'};    % S12 - Full  stim only experiment

% Set up paths
dataDir                = fullfile(fmsRootPath, 'data');    % Where to get data?
useSLIncohSpectrum     = true;      % Plot SL amplitudes from incoherent spectrum (default: true)

% Preallocate space for matrices
diffFullBlankSL = NaN(length(subject),157);
diffFullBlankBB = diffFullBlankSL;

%% 1. Load subject's data

for s = 1:length(subject)
    
    % Go from subject to session nr
    switch subject{s}
        case 'wlsubj002' % S1 - Full, Left, Right stim experiment
            whichSession = 2;
        case 'wlsubj004' % S2 - Full, Left, Right stim experiment
            whichSession = 7;
        case 'wlsubj005' % S3 - Full, Left, Right stim experiment
            whichSession = 8;
        case 'wlsubj006' % S4 - Full, Left, Right stim experiment
            whichSession = 1;
        case 'wlsubj010' % S5 - Full, Left, Right stim experiment
            whichSession = 6;
        case 'wlsubj011' % S6 - Full, Left, Right stim experiment
            whichSession = 5;
        case 'wlsubj048' % S7 - Full  stim only experiment
            whichSession = 9;
        case 'wlsubj046' % S8 - Full  stim only experiment
            whichSession = 10;
        case 'wlsubj039' % S9 - Full  stim only experiment
            whichSession = 11;
        case 'wlsubj059' % S10 - Full  stim only experiment
            whichSession = 12;
        case 'wlsubj067' % S11 - Full  stim only experiment
            whichSession = 13;
        case 'wlsubj070' % S12 - Full  stim only experiment
            whichSession = 14;
    end
    

    % Get amplitude data
    data = loadData(fullfile(dataDir, subject{s}), whichSession,'amplitudes');
    
    % Update SL amplitudes for each subject, either with coherent or incoherent spectrum
    if useSLIncohSpectrum
        ampl{s}.sl.full  = data.sl.full;
        ampl{s}.sl.blank = data.sl.blank;
        if strcmp(subject{s},'wlsubj059')
            ampl{s}.sl.full  = data.sl.full_coherent;
            ampl{s}.sl.blank = data.sl.blank_coherent;
        end
    else
        ampl{s}.sl.full  = data.sl.full_coherent;
        ampl{s}.sl.blank = data.sl.blank_coherent;
    end
    
    % Update broadband power for each subject, always from incoherent spectrum
    ampl{s}.bb.full  = data.bb.full;
    ampl{s}.bb.blank = data.bb.blank;
    
    clear data bb sl snr_sl snr_bb data;
    

    %% 2. Get contrast between full and blank for SL and BB data

    % Take difference between mean of full and blank epochs for each subject
    % and dataset (sl or bb)    
    diffFullBlankSL(s,:) = nanmean(ampl{s}.sl.full,1) - nanmean(ampl{s}.sl.blank,1);
    diffFullBlankBB(s,:) = nanmean(ampl{s}.bb.full,1) - nanmean(ampl{s}.bb.blank,1);
    
    % there is a difference between the datasets in terms of scaling units
    % session 1-6 are in fempto Tesla  whereas 7-12 are in Tesla
    % TODO: handle this more gracefully? Maybe just mark the sessions?
    if max(diffFullBlankSL(s,:)) < 1^-14
        diffFullBlankSL(s,:) = diffFullBlankSL(s,:) .* 10^15;
        diffFullBlankBB(s,:) = diffFullBlankBB(s,:) .* 10^15 .* 10^15;
    end

end




% template variables
eccenLimitDeg = [0 11];

% Define vector that can truncate number of sensors 
keep_sensors = logical([ones(157,1); zeros(192-157,1)]);

% Path to brainstorm database and project name
bsDB            = '/Volumes/server/Projects/MEG/brainstorm_db/';
projectName     = 'SSMEG';

for s = 1:length(subject)
    
    % Data point images
    
%     prediction = [diffFullBlankSL(s,:); diffFullBlankBB(s,:); zeros(1,157)];
%     prediction = [prediction, zeros(3,192-157)]';
%     prediction(:,2) = prediction(:,2).*10;
%     ok = savePrediction2SQD(fmsRootPath, subject{s}, prediction, [subject{s} '_dataSummary'])
    

    % Model point images
%     load(fullfile(dataDir, subject{s}, [subject{s} '_prediction_V123_0.00-11.mat'])) 
%     prediction = [dataToPlot zeros(3,192-157)]';
%     prediction = prediction.*10^4;
%     ok = savePrediction2SQD(fmsRootPath, subject{s}, prediction, [subject{s} '_modelSummary'])

    % Template V1-V3 within 11 deg eccentricity
    d = dir(fullfile(bsDB, projectName, 'data', subject{s}, 'R*'));
    bsData = fullfile(d(1).folder, d(1).name);    
    bsAnat = fullfile(bsDB, projectName, 'anat', subject{s});
    
    % Load Gain matrix
    G_constrained = getGainMatrix(bsData, keep_sensors);

    % Get V1 template limited to 11 degrees eccentricity
    template = getTemplate(bsAnat, 'V123', eccenLimitDeg);
    predictionV123 = template.V123_StimEccen*G_constrained';
    
    template = getTemplate(bsAnat, 'V1', eccenLimitDeg);
    predictionV1 = template.V1_StimEccen*G_constrained';

    template = getTemplate(bsAnat, 'V2', eccenLimitDeg);
    predictionV2 = template.V2_StimEccen*G_constrained';

    template = getTemplate(bsAnat, 'V3', eccenLimitDeg);
    predictionV3 = template.V3_StimEccen*G_constrained';
    
    prediction = [predictionV123', predictionV1', predictionV2', predictionV3'];
    prediction = [prediction; zeros(192-157,4)].*10^6;
    
    ok = savePrediction2SQD(fmsRootPath, subject{s}, prediction, [subject{s} '_modelV123_11deg'])
end



    
end
