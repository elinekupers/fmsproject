function makeSupplementaryFigure_coherentSpectrum11_to_13Hz(varargin)
% This is a function to make Supplementary Figure XX from the manuscript:
%   A visual encoding model links magnetoencephalography
%   signals to neural synchrony in human cortex.
%       by Kupers, Benson, Winawer (YEAR) JOURNAL.
%
% This figure shows the spatial pattern across MEG for 11, 12 and 13 Hz 
% for large field flickering (12 Hz) dartboard patterns, using the coherent
% spectrum.
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
%   [amplitudeType]     : (str)   plot amplitude at stimulus contrast 
%                                 reversal rate: 12 Hz (use 'amplitudes'),
%                                 or mean of amplitudes at the higher
%                                 frequencies used to define the broadband
%                                 power: 72,84,96,108,132,144 Hz (use
%                                 'amplitudesHigherHarmonics') (default is
%                                 'amplitudes')
%
% Example 1: Plot first subject
%  makeSupplementaryFigure_coherentSpectrum11_to_13Hz('subjectsToPlot', 1, 'plotMeanSubject', false, 'saveFig', true)
% Example 2: Plot example subject in manuscript (S12)
%  makeSupplementaryFigure_coherentSpectrum11_to_13Hz('subjectsToPlot', 12, 'plotMeanSubject', false, 'saveFig', true)
% Example 3: Plot all subjects and group average
%  makeSupplementaryFigure_coherentSpectrum11_to_13Hz('subjectsToPlot', 1:12, 'plotMeanSubject', true, 'saveFig', true)
%
% By Eline Kupers (NYU) 2017

%% 0. Set up paths and define parameters
p = inputParser;
p.KeepUnmatched = true;
p.addParameter('subjectsToPlot', 12);
p.addParameter('plotMeanSubject', true, @islogical); 
p.addParameter('saveFig', true, @islogical); 
p.addParameter('contourPercentile', [], @isnumeric);  
p.addParameter('maxColormapPercentile', 100, @isnumeric); 
p.addParameter('signedColorbar', false, @islogical);
p.addParameter('singleColorbarFlag', false, @islogical);
p.addParameter('amplitudeType', 'amplitudesCoherentSpectrum', ...
    @(x) any(validatestring(x,{'amplitudes', 'amplitudesHigherHarmonics', 'amplitudesCoherentSpectrum'})));
p.parse(varargin{:});

% Rename variables
subjectsToPlot          = p.Results.subjectsToPlot;
plotMeanSubject         = p.Results.plotMeanSubject;
saveFig                 = p.Results.saveFig;
contourPercentile       = p.Results.contourPercentile;
maxColormapPercentile   = p.Results.maxColormapPercentile;
signedColorbar          = p.Results.signedColorbar;
singleColorbarFlag      = p.Results.singleColorbarFlag;
amplitudeType           = p.Results.amplitudeType;

% Get subject names and corresponding data session number
[subject, dataSession] = getSubjectIDs;

% Set up paths
figureDir        = fullfile(fmsRootPath, 'figures_CoherentSpectrum'); % Where to save images?
dataDir          = fullfile(fmsRootPath, 'data');    % Where to get data?

% Plotting params
markerType    = '.';
colorContours = {'y','b'};

% Load all subjects when plotting the mean
if plotMeanSubject
    subjectsToLoad = 1:length(subject);
else
    subjectsToLoad = subjectsToPlot;
end

allData = cell(size(subjectsToLoad));

%% 1. Load subject's data
for s = subjectsToLoad
    
    % Go from subject to session nr
    whichSession = dataSession(s);
    
    % Get amplitude data
    data = loadData(fullfile(dataDir, subject{s}), whichSession, 'type', amplitudeType, 'plotNeighboringSLFreq', true);
    allData{s} = data;
    
    clear data
end

%% 3. Plot subject
for s = subjectsToPlot
    
    dataToPlotFull   = cat(1, allData{s}.sl.amps_full11, ...
                          allData{s}.sl.amps_full12, ...
                          allData{s}.sl.amps_full13);
    
    dataToPlotBlank   = cat(1, allData{s}.sl.amps_blank11, ...
                          allData{s}.sl.amps_blank12, ...
                          allData{s}.sl.amps_blank13);
               
    fig_ttlFull       = {sprintf('SupplFigureXX_Observed_MEGData_FullField_%s_prctile%2.1f_S%d_singleColorbarFlag%d', amplitudeType, contourPercentile, s, singleColorbarFlag), ...
                     sprintf('SupplFigureXX_Contour_FullField_%s_prctile%2.1f_S%d_singleColorbarFlag%d', amplitudeType, contourPercentile, s, singleColorbarFlag)};
    sub_ttl       = {sprintf('Amplitude at 11 Hz S%d', s), ...
                     sprintf('Amplitude at 12 Hz S%d', s), ...
                     sprintf('Amplitude at 13 Hz S%d', s)};
    
    if saveFig
        figureDirSubj = fullfile(figureDir, subject{s});
        if ~exist(figureDirSubj, 'dir'); mkdir(figureDirSubj); end
    else, figureDirSubj = []; 
    end
    
    % Plot full-field stimulus epochs
    visualizeSensormaps(dataToPlotFull, maxColormapPercentile, contourPercentile, ...
        signedColorbar, singleColorbarFlag, colorContours, markerType, fig_ttlFull, sub_ttl, saveFig, figureDirSubj);
    
    % Plot blank epochs
    fig_ttlBlank       = {sprintf('SupplFigureXX_Observed_MEGData_Blank_%s_prctile%2.1f_S%d_singleColorbarFlag%d', amplitudeType, contourPercentile, s, singleColorbarFlag), ...
                     sprintf('SupplFigureXX_Contour_Blank_%s_prctile%2.1f_S%d_singleColorbarFlag%d', amplitudeType, contourPercentile, s, singleColorbarFlag)};
                 
    visualizeSensormaps(dataToPlotBlank, maxColormapPercentile, contourPercentile, ...
        signedColorbar, singleColorbarFlag, colorContours, markerType, fig_ttlBlank, sub_ttl, saveFig, figureDirSubj);
   
end

%% 4. Plot average across subjects if requested

if plotMeanSubject
    
    for ii = subjectsToLoad
        amps11HzFull(1,ii,:) = allData{ii}.sl.amps_full11;
        amps12HzFull(1,ii,:) = allData{ii}.sl.amps_full12;
        amps13HzFull(1,ii,:) = allData{ii}.sl.amps_full13;
        
        amps11HzBlank(1,ii,:) = allData{ii}.sl.amps_blank11;
        amps12HzBlank(1,ii,:) = allData{ii}.sl.amps_blank12;
        amps13HzBlank(1,ii,:) = allData{ii}.sl.amps_blank13;
    end
    
    mnAmp.sl11Full = squeeze(mean(amps11HzFull(1,:,:),2,'omitnan'))';
    mnAmp.sl12Full = squeeze(mean(amps12HzFull(1,:,:),2,'omitnan'))';
    mnAmp.sl13Full = squeeze(mean(amps13HzFull(1,:,:),2,'omitnan'))';

    mnAmp.sl11Blank = squeeze(mean(amps11HzBlank(1,:,:),2,'omitnan'))';
    mnAmp.sl12Blank = squeeze(mean(amps12HzBlank(1,:,:),2,'omitnan'))';
    mnAmp.sl13Blank = squeeze(mean(amps13HzBlank(1,:,:),2,'omitnan'))';
    
    % Concatenate data
    dataToPlotFull      = cat(1, mnAmp.sl11Full, mnAmp.sl12Full, mnAmp.sl13Full);
    dataToPlotBlank      = cat(1, mnAmp.sl11Blank, mnAmp.sl12Blank, mnAmp.sl13Blank);

    
    % Define figure and subfigure titles
    fig_ttlFull         = {sprintf('SupplFigure11_to_13Hz_Observed_MEGData_FullField_%s_prctile%2.1f_singleColorbarFlag%d_AVERAGE', amplitudeType, contourPercentile,singleColorbarFlag), ...
                       sprintf('SupplFigure11_to_13Hz_Contour_%s_prctile%2.1f_singleColorbarFlag%d_AVERAGE', amplitudeType, contourPercentile,singleColorbarFlag)};
    sub_ttl         = {sprintf('Stimulus locked Average N = %d, 11 hz', length(subject)), ...
                     sprintf('Stimulus locked Average N = %d, 12 hz', length(subject)), ...
                     sprintf('Stimulus locked Average N = %d, 13 hz', length(subject))};
    
    % Make figure dir for average subject, if non-existing
    if saveFig
        figureDirAvg       = fullfile(figureDir,'average'); % Where to save images?
        if ~exist(figureDirAvg,'dir'); mkdir(figureDirAvg); end
    else, figureDirAvg = []; 
    end
    
    % Plot full field stimulus epochs
    visualizeSensormaps(dataToPlotFull, maxColormapPercentile, contourPercentile, ...
        signedColorbar, singleColorbarFlag, colorContours, markerType, fig_ttlFull, sub_ttl, saveFig, figureDirAvg);

    
     % Plot blank  epochs
     fig_ttlBlank         = {sprintf('SupplFigure11_to_13Hz_Observed_MEGData_Blank_%s_prctile%2.1f_singleColorbarFlag%d_AVERAGE', amplitudeType, contourPercentile,singleColorbarFlag), ...
                       sprintf('SupplFigure11_to_13Hz_Contour_%s_prctile%2.1f_singleColorbarFlag%d_AVERAGE', amplitudeType, contourPercentile,singleColorbarFlag)};
    
    % Plot full field stimulus epochs
    visualizeSensormaps(dataToPlotBlank, maxColormapPercentile, contourPercentile, ...
        signedColorbar, singleColorbarFlag, colorContours, markerType, fig_ttlBlank, sub_ttl, saveFig, figureDirAvg);

     
end

return



