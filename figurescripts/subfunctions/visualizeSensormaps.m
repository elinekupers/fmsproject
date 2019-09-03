function sensorsOfInterest = visualizeSensormaps(data, maxColormapPercentile, contourPercentile, signedColorbar, colorMarkers, markerType, fig_ttl, sub_ttl, saveFigures, figureDir)

% visualizeSensormaps(data, colormapLims, contourmapLims, colorMarkers, markerType, fig_ttl, sub_ttl, saveFigures)

% Function to visualize 4 sensormaps and the contour lines based on the data
% of the 2 rows

% Example:
% data = rand(4,157);
% visualizeSensormaps(data, 50, 50)

%% Check input params
if ~exist('maxColormapPercentile','var') || isempty(maxColormapPercentile)
    maxColormapPercentile = 100;
end

if ~exist('contourPercentile','var') || isempty(contourPercentile)
    contourPercentile = []; % if empty, don't show any
end

if ~exist('signedColorbar','var') || isempty(signedColorbar)
    signedColorbar = true; 
end

if ~exist('colorMarkers','var') || isempty(colorMarkers)
    if size(data,1) == 3
        colorMarkers =  repmat({'y','b','g'},1,2);
    else
        colorMarkers =  repmat({'y','b'},1,2);
    end
end

if ~exist('markerType','var') || isempty(markerType)
    markerType = '.'; % could also be '*';
end

if ~exist('fig_ttl','var') || isempty(fig_ttl)
    fig_ttl = {'Figure1','Figure2'};
end

if ~exist('sub_ttl','var') || isempty(sub_ttl)
     sub_ttl = {'','',''};
end

if ~exist('saveFigures','var') || isempty(saveFigures)
    saveFigures = false;
    figureDir = [];
else
    if saveFigures && (~exist('figureDir','var') || isempty(saveFigures))
        figureDir = fullfile(fmsRootPath, 'figures'); 
        if ~exist(figureDir,'dir'); mkdir(figureDir); end
    end
end

%% Predefine figures and matrix saving sensor nrs drawn within contours

% Main figure
fH1 = figure; clf; set(fH1,'position',[1,600,1400,800], 'Name', fig_ttl{1}, 'NumberTitle', 'off');

% If contour lines are getting drawn, replot them in second figure
if ~isempty(contourPercentile) && contourPercentile>10
    cmapContour = [1 1 1; 0 0 1; 1 0 0];
    fH2 = figure; clf; set(fH2,'position',[ 1400, 923, 700, 473], 'Name', fig_ttl{2}, 'NumberTitle', 'off');
                   subplot(1,1,1); megPlotMap(zeros(1,157)); colormap(cmapContour);
end

% Save sensors of interest that fall within the contour lines
sensorsOfInterest = NaN(size(data));

%% Loop over datasets
for ii = 1:size(data,1)
    
    dataToPlot = data(ii,:);
    cmapData   = bipolar(64);
    
    % Check limits of color map/bar 
    if signedColorbar
        colormapLims = [-1,1].*prctile(dataToPlot, maxColormapPercentile);
    else
        colormapLims = [0 prctile(dataToPlot, maxColormapPercentile)];
        cmapData = cmapData(ceil(length(cmapData)/2):end,:);
    end

    % Check at what point to draw contour lines
    if ~isempty(contourPercentile)
       if contourPercentile > 10 % Plot at given data percentile 
           contourLines = [1 1]*prctile(dataToPlot, contourPercentile); 
       else % Plot nr of contours at equo-distance percentiles 
           contourLines = contourPercentile; 
       end
    else, contourLines = []; end


    % Plot data or predictions
    figure(fH1);
    subplot(size(data,1),1,ii);
    [~,ch] = megPlotMap(dataToPlot,colormapLims,fH1,cmapData,sub_ttl{ii},[],[], ...
        'isolines', contourLines, ...
    ...    'chanindx', dataToPlot > max(contourmapLims), ...
        'pointsymbol', markerType, ... '*'
        'pointsize', 10); hold on;
    % Check if a contour line was requested, if so, save those data and
    % replot in separate figure
    if ~isempty(contourLines)
        c = findobj(gca,'Type','Contour');
%         pp = findobj(gca,'Marker',markerType);
        
        for jj = 1:size(c,1)
            % Change line width
            figure(fH1); 
            c(jj).LineWidth = 2;
        end
        
        if length(contourLines)==2
            figure(fH2); hold all;
            
            contourf(c.XData, c.YData, c.ZData, contourLines, 'LineColor',colorMarkers{ii}, 'Fill','off','LineWidth',4);
            %%scatter(pp(1).XData,pp(1).YData, 150, colorMarkers{ii},'*'); 
            colorbar off;
            
            sensorsOfInterest(ii,:) = dataToPlot > max(contourLines(1));
        end
    end
    
    % Make figure pretty
    set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12); title(sub_ttl{ii})
    
end

% Save if requested
if saveFigures
    if ~exist('figureDir', 'dir'); mkdir(figureDir); end
    set(0, 'currentfigure', fH1);
    figurewrite(fullfile(figureDir,fig_ttl{1}),[],0,'.',1);
    if length(contourLines)==2
        set(0, 'currentfigure', fH2);
        figurewrite(fullfile(figureDir,fig_ttl{2}),[],0,'.',1);
    end
end