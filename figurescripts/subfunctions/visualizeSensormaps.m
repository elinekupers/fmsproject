function sensorsOfInterest = visualizeSensormaps(data, colormapPercentile, contourmapPercentile, colorMarkers, markerType, fig_ttl, sub_ttl, saveFigures, figureDir)

% visualizeSensormaps(data, colormapLims, contourmapLims, colorMarkers, markerType, fig_ttl, sub_ttl, saveFigures)

% Function to visualize 4 sensormaps and the contour lines based on the data
% of the 2 rows

% Example:
% data = rand(4,157);
% visualizeSensormaps(data, 50, 50)

%% Check input params
if ~exist('colormapPercentile','var') || isempty(colormapPercentile)
    colormapPercentile = 97.5;
end

if ~exist('contourmapPercentile','var') || isempty(contourmapPercentile)
    contourmapPercentile = 93.6; % such that we show the top 10 sensors
end

if ~exist('colorMarkers','var') || isempty(colorMarkers)
    colorMarkers = {'r','b', 'r', 'b'};
end

if ~exist('markerType','var') || isempty(markerType)
    markerType = '.'; % could also be '*';
end

if ~exist('fig_ttl','var') || isempty(fig_ttl)
    fig_ttl = {'Figure1','Figure2'};
end

if ~exist('sub_ttl','var') || isempty(sub_ttl)
    sub_ttl = {'','','',''};
end

if ~exist('saveFigures','var') || isempty(saveFigures)
    saveFigures = false;
    figureDir = [];
else
    if saveFigures && (~exist('figureDir','var') || isempty(saveFigures))
        figureDir = fullfile(fmsRootPath, 'figures'); if ~exist(figureDir,'dir'); mkdir(figureDir); end
    end
end

%% Predefine figures 
fH1 = figure; clf; set(fH1,'position',[1,600,1400,800], 'Name', fig_ttl{1}, 'NumberTitle', 'off');
fH2 = figure; clf; set(fH2,'position',[1400,600,700,800], 'Name', fig_ttl{2}, 'NumberTitle', 'off');
                   subplot(2,1,1); megPlotMap(zeros(1,157)); colormap([1 1 1]);
                   subplot(2,1,2); megPlotMap(zeros(1,157)); colormap([1 1 1]);

% Save sensors of interest that fall within the contour lines
sensorsOfInterest = NaN(size(data));

% Loop over datasets
for ii = 1:size(data,1)
    
    dataToPlot = data(ii,:);
    colormapLims =  [-1 1]*prctile(dataToPlot, colormapPercentile);
    contourmapLims = [1 1]*prctile(dataToPlot, contourmapPercentile);

    % Plot data or predictions
    figure(fH1);
    subplot(2,2,ii);
    [~,ch] = megPlotMap(dataToPlot,colormapLims,fH1,'bipolar',[],[],[], ...
        'isolines', contourmapLims, ...
    ...    'chanindx', dataToPlot > max(contourmapLims), ...
        'pointsymbol', markerType, ... '*'
        'pointsize', 10);
    
    c = findobj(gca,'Type','Contour'); c.LineWidth = 4;
    pp = findobj(gca,'Marker',markerType);
    set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12); title(sub_ttl{ii})
    
    % Plot overlap
    figure(fH2);
    subplot(2,1,ceil(ii/2)); hold all;   
    contour(c.XData, c.YData, c.ZData, contourmapLims, 'LineColor',colorMarkers{ii}, 'LineWidth',4);
    %%scatter(pp(1).XData,pp(1).YData, 150, colorMarkers{ii},'*'); 
    colorbar off;
    
    sensorsOfInterest(ii,:) = dataToPlot > max(contourmapLims);
    
end

% Save if requested
if saveFigures
    set(0, 'currentfigure', fH1);
    figurewrite(fullfile(figureDir,fig_ttl{1}),[],0,'.',1);
    set(0, 'currentfigure', fH2);
    figurewrite(fullfile(figureDir,fig_ttl{2}),[],0,'.',1);
end