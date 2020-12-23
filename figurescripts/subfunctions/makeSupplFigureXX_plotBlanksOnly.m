function makeSupplFigureXX_plotBlanksOnly(varargin)
% This is a function to make Supplementary Figure XX from the manuscript:
%   A visual encoding model links magnetoencephalography
%   signals to neural synchrony in human cortex.
%       by Kupers, Benson, Winawer (YEAR) JOURNAL.
%
% This figure shows the spatial pattern across MEG for the blank periods 
% shown in between the large field flickering (12 Hz) dartboard pattern.
%
% To runs this script, you need:
% (1) the data from the denoiseproject in the data folder of its FMS
%     code repository
%
% (2) MEG_utils toolbox added to the paths. For example:
%     tbUse('ForwardModelSynchrony');
%        or to only add the MEG_utils toolbox:
%     addpath(genpath('~/matlab/git/toolboxes/meg_utils'))
%
% INPUTS:
%   [subjectsToPlot]     : (int)  subject nr to plot (default: 12)
%   [plotMeanSubject]    : (bool) plot average across all 12 subjets or not?
%                                 (default: true)
%   [saveFig]            : (bool) save figures or not? (default: true)
%   [useSLIncohSpectrum] : (bool) plot SL amplitudes from incoherent spectrum
%                                 (default: true)
%   [contourPercentile]  : (int)  percentile of the data to draw contour 
%                                 lines? There are two ways to define: 
%                                 (1) single integer above 10 to get percentile
%                                 to select X sensors with highest response
%                                 Use 90.4 for top 15 sensors, or 93.6 for
%                                 top 10 sensors.
%                                 (2) single integer below 10 to get contour 
%                                 lines at equal percentiles of data, e.g.
%                                 3 => 25 50 75th prctle
%                                 (default: 93.6) 
%   [maxColormapPercentile]:(int) percentile of data to truncate colormap
%                                 (default: 97.5)
%   [signedColorbar]    : (bool)  plot signed colormap (true) or only 
%                                 positive values (false)? 
%                                 (default: true)
%   [singleColorbarFlag]: (bool)  Use a single colorbar per row,instead of 
%                                 per individual plot (default: false)
%   [snrThresh]         : (int)   snr threshold for amplitudes
%                                 (default: 1)
%   [useSLPower]        : (bool)  plot SL power or amplitudes from spectrum
%                                 (default: false)
%   [amplitudeType]     : (str)   plot amplitude at stimulus contrast 
%                                 reversal rate: 12 Hz (use 'amplitudes'),
%                                 or mean of amplitudes at the higher
%                                 frequencies used to define the broadband
%                                 power: 72,84,96,108,132,144 Hz (use
%                                 'amplitudesHigherHarmonics') (default is
%                                 'amplitudes')
%   [fixedRNGseed]      : (bool)  use a fixed random number generator seed,
%                                 for computational reproducibility.
%                                 (default: false)
%
% Example 1: Plot first subject
%  makeSupplFigureXX_plotBlanksOnly('subjectsToPlot', 1, 'plotMeanSubject', false, 'saveFig', true)
% Example 2: Plot example subject in manuscript (S12)
%  makeSupplFigureXX_plotBlanksOnly('subjectsToPlot', 12, 'plotMeanSubject', false, 'saveFig', true)
% Example 3: Plot all subjects and group average
%  makeSupplFigureXX_plotBlanksOnly('subjectsToPlot', 1:12, 'plotMeanSubject', true, 'saveFig', true)
% Example 4: Plot all subjects and group average with fixed rng seed
%  makeSupplFigureXX_plotBlanksOnly('subjectsToPlot', 1:12, 'plotMeanSubject', true, 'saveFig', true, 'fixedRNGseed', true)
%
% By Eline Kupers (NYU) 2017

%% 0. Set up paths and define parameters
p = inputParser;
p.KeepUnmatched = true;
p.addParameter('subjectsToPlot', 12);
p.addParameter('plotMeanSubject', true, @islogical); 
p.addParameter('saveFig', true, @islogical); 
p.addParameter('useSLIncohSpectrum', true, @islogical);
p.addParameter('contourPercentile', [], @isnumeric);  
p.addParameter('maxColormapPercentile', 97.5, @isnumeric); 
p.addParameter('signedColorbar', false, @islogical);
p.addParameter('singleColorbarFlag', false, @islogical);
p.addParameter('snrThresh',1, @isnumeric);
p.addParameter('useSLPower', false, @islogical);
p.addParameter('amplitudeType', 'amplitudes', ...
    @(x) any(validatestring(x,{'amplitudes', 'amplitudesHigherHarmonics'})));
p.addParameter('fixedRNGseed', false, @islogical);
p.parse(varargin{:});

% Rename variables
subjectsToPlot          = p.Results.subjectsToPlot;
plotMeanSubject         = p.Results.plotMeanSubject;
saveFig                 = p.Results.saveFig;
useSLIncohSpectrum      = p.Results.useSLIncohSpectrum;
contourPercentile       = p.Results.contourPercentile;
maxColormapPercentile   = p.Results.maxColormapPercentile;
signedColorbar          = p.Results.signedColorbar;
singleColorbarFlag      = p.Results.singleColorbarFlag;
snrThresh               = p.Results.snrThresh;
useSLPower              = p.Results.useSLPower;
amplitudeType           = p.Results.amplitudeType;
fixedRNGseed            = p.Results.fixedRNGseed;

% if requested: set a fixed random number generator seed (for
% computational reproducibility of selecting a random subset of epochs 
% for some subjects in loadData.m).
if fixedRNGseed
    rng(1);
end

% Get subject names and corresponding data session number
[subject, dataSession] = getSubjectIDs;

% Set up paths
figureDir        = fullfile(fmsRootPath, 'figuresBlankOnly'); % Where to save images?
dataDir          = fullfile(fmsRootPath, 'data');    % Where to get data?

% Plotting params
markerType    = '.';
colorContours = {'y','b'};

% Preallocate space for matrices
fullBlankSL = NaN(length(subject),157);
fullBlankBB = fullBlankSL;

% Load all subjects when plotting the mean
if plotMeanSubject
    subjectsToLoad = 1:length(subject);
else
    subjectsToLoad = subjectsToPlot;
end

ampl = cell(size(subjectsToLoad));

%% 1. Load subject's data
for s = subjectsToLoad
    
    % Go from subject to session nr
    whichSession = dataSession(s);
    
    % Get SNR data
    data = loadData(fullfile(dataDir, subject{s}),whichSession,'SNR');
    SNR(s,1,:) = data{1}; %#ok<AGROW>
    SNR(s,2,:) = data{2}; %#ok<AGROW>
    
    % Get amplitude data
    data = loadData(fullfile(dataDir, subject{s}), whichSession, amplitudeType, fixedRNGseed);
    
    % Update SL amplitudes for each subject, either with coherent or incoherent spectrum
    if useSLIncohSpectrum 
        ampl{s}.sl.blank = data.sl.blank;
        if (~strcmp(amplitudeType, 'amplitudesHigherHarmonics') && strcmp(subject{s},'wlsubj059'))
            ampl{s}.sl.blank = data.sl.blank_coherent;
        end
    else
        ampl{s}.sl.blank = data.sl.blank_coherent;
    end
    
    % Convert SL amplitudes to power (by squaring) if requested
    if useSLPower
        ampl{s}.sl.blank = ampl{s}.sl.blank.^2;
    end
    
    % Update broadband power for each subject, always from incoherent spectrum;
    ampl{s}.bb.blank = data.bb.blank;
    
    clear data bb sl snr_sl snr_bb data;

    %% 2. Get contrast between full and blank for SL and BB data

    % Take difference between mean of full and blank epochs for each subject
    % and dataset (sl or bb)
    fullBlankSL(s,:) = nanmean(ampl{s}.sl.blank,1);
    fullBlankBB(s,:) = nanmean(ampl{s}.bb.blank,1);
end

for s = subjectsToPlot
    %% 3. Plot subject
    snrThreshMask.sl.single = abs(squeeze(SNR(s,1,:))) > snrThresh;
    snrThreshMask.bb.single = abs(squeeze(SNR(s,2,:))) > snrThresh;

    dataToPlot   = cat(1, fullBlankSL(s,:) .* snrThreshMask.sl.single', ...
       fullBlankBB(s,:) .* snrThreshMask.bb.single');
     
    fig_ttl       = {sprintf('Figure3_Observed_MEG_BlankData_incohSpectrum%d_prctile%2.1f_S%d_slPower%d_singleColorbarFlag%d', useSLIncohSpectrum, contourPercentile, s, useSLPower,singleColorbarFlag), ...
                     sprintf('Figure3_Contour_incohSpectrum%d_prctile%2.1f_S%d_slPower%d_singleColorbarFlag%d', useSLIncohSpectrum, contourPercentile, s, useSLPower,singleColorbarFlag)};
    sub_ttl       = {sprintf('Stimulus locked S%d', s), ...
                     sprintf('Broadband S%d', s)};
    
    if saveFig
        figureDirSubj = fullfile(figureDir, subject{s});
        if ~exist(figureDirSubj, 'dir'); mkdir(figureDirSubj); end
    else, figureDirSubj = []; 
    end
    
    % Plot it!
    visualizeSensormaps(dataToPlot, maxColormapPercentile, contourPercentile, ...
        signedColorbar, singleColorbarFlag, colorContours, markerType, fig_ttl, sub_ttl, saveFig, figureDirSubj);
end

%% 4. Plot average across subjects if requested

if plotMeanSubject
    
    % Get snr mask for group data
    snrThreshMask.sl.group  = abs(squeeze(nanmean(SNR(:,1,:),1))) > snrThresh;
    snrThreshMask.bb.group  = abs(squeeze(nanmean(SNR(:,2,:),1))) > snrThresh;
    
    % Concatenate data
    dataToPlot      = cat(1, nanmean(fullBlankSL,1) .* snrThreshMask.sl.group', ...
                             nanmean(fullBlankBB,1) .* snrThreshMask.bb.group');
    
    % Define figure and subfigure titles
    fig_ttl         = {sprintf('Figure3_Observed_MEG_BlankData_incohSpectrum%d_prctile%2.1f_slPower%d_singleColorbarFlag%d_AVERAGE', useSLIncohSpectrum, contourPercentile, useSLPower,singleColorbarFlag), ...
                       sprintf('Figure3_Contour_incohSpectrum%d_prctile%2.1f_slPower%d_singleColorbarFlag%d_AVERAGE', useSLIncohSpectrum, contourPercentile, useSLPower,singleColorbarFlag)};
    sub_ttl         = {sprintf('Stimulus locked Average N = %d', length(subject)), ...
                       sprintf('Broadband Average N = %d', length(subject))};
    
    % Make figure dir for average subject, if non-existing
    if saveFig
        figureDirAvg       = fullfile(figureDir,'average'); % Where to save images?
        if ~exist(figureDirAvg,'dir'); mkdir(figureDirAvg); end
    else, figureDirAvg = []; 
    end
    
    % Plot it!
    visualizeSensormaps(dataToPlot, maxColormapPercentile, contourPercentile, ...
        signedColorbar, singleColorbarFlag, colorContours, markerType, fig_ttl, sub_ttl, saveFig, figureDirAvg);
end

return



