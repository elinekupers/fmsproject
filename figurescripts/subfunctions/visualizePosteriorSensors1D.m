function [] = visualizePosteriorSensors1D(data, plotAvg, fig_ttl,sub_ttl, saveFigures, figureDir)
% Function to visualize 1D summary of posterior MEG sensors data 
% from the manuscript:
%   A visual encoding model links magnetoencephalography
%   signals to neural synchrony in human cortex.
%       by Kupers, Benson, Winawer (YEAR) JOURNAL.
%
% This figure shows the weighted sum of MEG sensor responses, using a 
% Gaussian pooling window in small bins going from left to right of the head.
% This gives us a 1D representation of the spatial pattern in the back half
% of the subjects head.
%
% INPUTS:
%   data                 : (struct or cell) sl and bb struct from
%                                 loadData.m
%   plotAvg              : (bool) plot average across all 12 subjets or not?
%                                 (default: false)
%   [fig_ttl]            : (str)  Figure title, used when saving figure.
%   [sub_ttl]            : (cell) Two subplot titles default: {'SL','BB'}
%   [saveFigures]        : (bool) save figures or not? (default: false)
%   [figureDir]          : (str)  folder to save figures
%
% Example 1:
% [subject, dataSession] = getSubjectIDs;
% subj = 1;
% data = loadData(fullfile(fmsRootPath, 'data', subject{subj}), dataSession(subj), 'type', 'amplitudes');
% visualizePosteriorSensors1D(data, false)
%
%% Handle inputs
if ~exist('plotAvg','var') || isempty(plotAvg)
    plotAvg = false;
end

if ~exist('fig_ttl','var') || isempty(fig_ttl)
    fig_ttl = {'1D MEG sensor data'};
end

if ~exist('sub_ttl','var') || isempty(sub_ttl)
    sub_ttl = {'SL', 'BB'};
end

if ~exist('saveFigures','var') || isempty(saveFigures)
    saveFigures = false;
end

if ~exist('figureDir','var') || isempty(figureDir)
    figureDir = [];
end

% get sensor locations
load(which('meg160_example_hdr.mat'))
layout = ft_prepare_layout([],hdr);
xpos   = layout.pos(1:157,1);
ypos   = layout.pos(1:157,2);

% Create x-pos bins and define width of gaussian summing window
nbins    = 100;
xbins    = linspace(min(xpos),max(xpos),nbins);
windowsz = 0.1;
sigma    = 0.05;

% plotting params
colors = {'r', 'b'};

% Put single subject data in a cell.
if ~iscell(data)
    tmp = data; clear data;
    data{1} = tmp;
end
% Get dimensions
[nBoot, ~, ~] = size(data{1}.sl.boot_amps_diff_mn);
nSubjects     = length(data);
nDataTypes    = 2;
    
% Preallocate space
mnBoot     = NaN(nSubjects,nBoot,nDataTypes,size(xbins,2));
sampleMean = NaN(nSubjects,nDataTypes,size(xbins,2));

for s = 1:nSubjects
    
    bootData = cat(3, data{s}.sl.boot_amps_diff_mn, data{s}.bb.boot_amps_diff_mn);
    
    for boot = 1:nBoot
        for dt = 1:nDataTypes
            thisSubjectData = squeeze(bootData(boot,:,dt));
            for c = 1:nbins
                %                 keep = abs(xpos-centers(c)) < 0.5*windowsz & ypos < 0;
                wt = exp(-0.5*((xbins(c)-xpos)/sigma).^2);
                wt(ypos>0)=0;
                mnBoot(s,boot,dt,c) = sum(thisSubjectData.*wt', 'omitnan')/sum(wt);
            end
        end
    end
    
    err = std(mnBoot, [],2,'omitnan');
    
    for c = 1:nbins
        %         keep = abs(xpos-centers(c)) < 0.5*windowsz & ypos < 0;
        wt = exp(-0.5*((xbins(c)-xpos)/sigma).^2);
        wt(ypos>0)=0;
        sampleMean(s,1,c) = sum(data{s}.sl.amps_diff_mn.*wt', 'omitnan')/sum(wt);
        sampleMean(s,2,c) = sum(data{s}.bb.amps_diff_mn.*wt', 'omitnan')/sum(wt);
    end
end

% Visualize data
xlbl = xbins;

if ~plotAvg
    figure; cla; set(gcf, 'color', 'w');
    if nSubjects > 6; nrows = 4; ncols = 6; else, nrows = 2; ncols = nSubjects; end
    
    for s = 1:nSubjects
        for ii = 1:nDataTypes
            subplot(nrows,ncols,s+(s*(ii-1)));
            plot(xlbl,squeeze(sampleMean(s,ii,:)), 'color', colors{ii},'LineWidth',6); hold on;
            errorbar(xlbl,squeeze(sampleMean(s,ii,:)), squeeze(err(s,1,ii,:)), 'k')
            if any(squeeze((sampleMean(s,ii,:)-min(err(s,1,ii,:))))<0)
                plot(xlbl,zeros(size(xlbl)), 'k','LineWidth',1);
                yl = [min(sampleMean(s,ii,:))-max(err(s,1,ii,:)), max(sampleMean(s,ii,:))+max(err(s,1,ii,:))];
            else
                yl = [0, max(sampleMean(s,ii,:))+max(err(s,1,ii,:))];
            end
            box off; title(sub_ttl{ii})
            set(gca, 'TickDir','out','ticklength',[0.010 0.010], 'FontSize',15);
            ylabel('Amplitude (fT)')
            ylim(yl);
        end
    end
    
    % Save if requested
    if saveFigures
        figurewrite(fullfile(figureDir,[fig_ttl '_individual']),[],0,'.',1);
    end
    
else
    
    groupMn = squeeze(mean(sampleMean,1,'omitnan'));
    groupErr = squeeze(std(sampleMean,[],1, 'omitnan')/sqrt(nSubjects));
    
    figure; cla; set(gcf, 'color', 'w');
    for ii = 1:nDataTypes
        subplot(2,1,ii);
        plot(xlbl,groupMn(ii,:), '-','color',colors{ii},'LineWidth',6); hold on;
        errorbar(xlbl,groupMn(ii,:), groupErr(ii,:), 'k')
        if any((groupMn(ii,:)-min(groupErr(ii,:)))<0)
            plot(xlbl,zeros(size(xlbl)), 'k','LineWidth',1);
            yl = [min(groupMn(ii,:))-min(groupErr(ii,:)), max(groupMn(ii,:))+max(groupErr(ii,:))];
        else
            yl = [0, max(groupMn(ii,:))+max(groupErr(ii,:))];
        end
        box off; title(sprintf('S%d', s))
        set(gca, 'TickDir','out','ticklength',[0.010 0.010], 'FontSize',15);
        ylabel('Amplitude (fT)')
        ylim(yl);
    end
    
    % Save if requested
    if saveFigures
        figurewrite(fullfile(figureDir,[fig_ttl '_GroupAverage']),[],0,'.',1);
    end
end



end