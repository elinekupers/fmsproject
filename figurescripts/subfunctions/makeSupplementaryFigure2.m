function makeSupplementaryFigure2(varargin)
% This is a function to make Supplementary Figure 2 from the manuscript:
%   A visual encoding model links magnetoencephalography
%   signals to neural synchrony in human cortex.
%       by Kupers, Benson, Winawer (YEAR) JOURNAL.
%
% This figure shows the asynchronous broadband response to a large field
% flickering (12 Hz) dartboard pattern, summarized as the geometric mean
% of 10 Hz bins between 60-150 Hz, excl. stimulus-locked harmonics
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
%
% Example 1: Plot first subject
%  makeSupplementaryFigure2('subjectsToPlot', 1, 'plotMeanSubject', false, 'saveFig', true)
% Example 2: Plot example subject in manuscript (S12)
%  makeSupplementaryFigure2('subjectsToPlot', 12, 'plotMeanSubject', false, 'saveFig', true)
% Example 3: Plot all subjects and group average
%  makeSupplementaryFigure2('subjectsToPlot', 1:12, 'plotMeanSubject', true, 'saveFig', true)
%
% By Eline Kupers (NYU) 2017

%% 0. Set up paths and define parameters
p = inputParser;
p.KeepUnmatched = true;
p.addParameter('subjectsToPlot', 12);
p.addParameter('plotMeanSubject', true, @islogical);
p.addParameter('saveFig', true, @islogical);
p.addParameter('contourPercentile', 93.6, @isnumeric);
p.addParameter('maxColormapPercentile', 97.5, @isnumeric);
p.addParameter('signedColorbar', true, @islogical);
p.addParameter('singleColorbarFlag', false, @islogical);
p.addParameter('snrThresh',0.5, @isnumeric);
p.addParameter('amplitudeType', 'amplitudesBB10hz');
p.parse(varargin{:});

% Rename variables
subjectsToPlot          = p.Results.subjectsToPlot;
plotMeanSubject         = p.Results.plotMeanSubject;
saveFig                 = p.Results.saveFig;
contourPercentile       = p.Results.contourPercentile;
maxColormapPercentile   = p.Results.maxColormapPercentile;
signedColorbar          = p.Results.signedColorbar;
singleColorbarFlag      = p.Results.singleColorbarFlag;
snrThresh               = p.Results.snrThresh;
amplitudeType           = p.Results.amplitudeType;

% Get subject names and corresponding data session number
[subject, dataSession] = getSubjectIDs;

% Set up paths
figureDir        = fullfile(fmsRootPath, 'figures'); % Where to save images?
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
    data = loadData(fullfile(dataDir, subject{s}), whichSession, 'type', amplitudeType);
    allData{s} = data;
    
    clear data
    
end

%% 3. Mask data by SNR, plot data per subject
figureDirs = {'_BB60-70Hz','_BB70-80Hz','_BB80-90Hz','_BB90-100Hz','_BB100-110Hz', ...
    '_BB110-120Hz','_BB120-130Hz','_BB130-140Hz','_BB140-150Hz'};


for s = subjectsToPlot
    for ii = 1:size(allData{s}.bb.amps_full,3)
        snrThreshMask.sl.single = abs(allData{s}.sl.snr) > snrThresh;
        snrThreshMask.bb.single = abs(allData{s}.bb.snr(:,:,ii)) > snrThresh;
        
        dataToPlot   = cat(1, allData{s}.sl.amps_diff_mn .* snrThreshMask.sl.single, ...
            allData{s}.bb.amps_diff_mn(:,:,ii) .* snrThreshMask.bb.single);
        
        fig_ttl       = {sprintf('SupplFig2_Observed_MEG_Data_%s_prctile%2.1f_S%d_singleColorbarFlag%d', amplitudeType, contourPercentile, s, singleColorbarFlag), ...
            sprintf('SupplFig2_Contour_%s_prctile%2.1f_S%d_singleColorbarFlag%d', amplitudeType, contourPercentile, s, singleColorbarFlag)};
        sub_ttl       = {sprintf('Stimulus locked S%d', s), ...
            sprintf('Broadband S%d', s)};
        
        if saveFig
            figureDirSubj = fullfile(figureDir,['figure' figureDirs{ii}], subject{s});
            if ~exist(figureDirSubj, 'dir'); mkdir(figureDirSubj); end
        else, figureDirSubj = [];
        end
        
        % Plot it!
        visualizeSensormaps(dataToPlot, maxColormapPercentile, contourPercentile, ...
            signedColorbar, singleColorbarFlag, colorContours, markerType, fig_ttl, sub_ttl, saveFig, figureDirSubj);
    end
    
    %     for s = subjectsToPlot
    %         figureDirSubj = fullfile(figureDir,['figure' figureDirs{ii}], subject{s});
    %         sub_ttl       = {sprintf('Stimulus locked S%d', s), ...
    %             sprintf('Broadband S%d', s)};
    %
    %         fig_ttl2 = sprintf('Figure3_1Daverage_Observed_MEG_Data_%s_S%d', amplitudeType, s);
    %         visualizePosteriorSensors1D(allData{s}, false, fig_ttl2,sub_ttl, saveFig, figureDirSubj)
    %     end
    
end

%% 4. Plot average across subjects if requested

if plotMeanSubject
    for jj = 1:size(allData{s}.bb.amps_full,3)
        for ii = subjectsToLoad
            snr(1,ii,:) = allData{ii}.sl.snr;
            snr(2,ii,:) = allData{ii}.bb.snr(:,:,jj);
            
            amps(1,ii,:) = allData{ii}.sl.amps_diff_mn;
            amps(2,ii,:) = allData{ii}.bb.amps_diff_mn(:,:,jj);
        end
        
        mnSNR.sl = squeeze(mean(snr(1,:,:),2,'omitnan'))';
        mnSNR.bb = squeeze(mean(snr(2,:,:),2,'omitnan'))';
        
        mnAmp.sl = squeeze(mean(amps(1,:,:),2,'omitnan'))';
        mnAmp.bb = squeeze(mean(amps(2,:,:),2,'omitnan'))';
        
        % Get snr mask for group data
        snrThreshMask.sl.group  = abs(mnSNR.sl) > snrThresh;
        snrThreshMask.bb.group  = abs(mnSNR.bb) > snrThresh;
        
        % Concatenate data
        dataToPlot      = cat(1, mnAmp.sl .* snrThreshMask.sl.group, ...
            mnAmp.bb .* snrThreshMask.bb.group);
        
        % Define figure and subfigure titles
        fig_ttl         = {sprintf('SupplFig2_Observed_MEG_Data_%s_prctile%2.1f_singleColorbarFlag%d_AVERAGE', amplitudeType, contourPercentile,singleColorbarFlag), ...
            sprintf('SupplFig2_Contour_%s_prctile%2.1f_singleColorbarFlag%d_AVERAGE', amplitudeType, contourPercentile,singleColorbarFlag)};
        sub_ttl         = {sprintf('Stimulus locked Average N = %d', length(subject)), ...
            sprintf('Broadband Average N = %d', length(subject))};
        
        % Make figure dir for average subject, if non-existing
        if saveFig
            figureDirAvg       = fullfile(figureDir,['figure' figureDirs{jj}],'average'); % Where to save images?
            if ~exist(figureDirAvg,'dir'); mkdir(figureDirAvg); end
        else, figureDirAvg = [];
        end
        
        % Plot it!
        visualizeSensormaps(dataToPlot, maxColormapPercentile, contourPercentile, ...
            signedColorbar, singleColorbarFlag, colorContours, markerType, fig_ttl, sub_ttl, saveFig, figureDirAvg);
    end
    
    fig_ttl2 = sprintf('SupplFig2_1Daverage_Observed_MEG_Data_%s_AVERAGE', amplitudeType);
    visualizePosteriorSensors1D(allData, true, fig_ttl2, sub_ttl,saveFig, figureDirAvg)
end

return



