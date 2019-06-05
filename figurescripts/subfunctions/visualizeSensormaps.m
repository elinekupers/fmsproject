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
    contourmapPercentile = []; % if empty, don't show any
end

if ~exist('colorMarkers','var') || isempty(colorMarkers)
    colorMarkers =  repmat({'r','b','g'},1,2);
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
cmap = [1 1 1; 0 0 1; 1 0 0];
fH1 = figure; clf; set(fH1,'position',[1,600,1400,800], 'Name', fig_ttl{1}, 'NumberTitle', 'off');

if ~isempty(contourmapPercentile) && contourmapPercentile>10
    fH2 = figure; clf; set(fH2,'position',[ 1400, 923, 700, 473], 'Name', fig_ttl{2}, 'NumberTitle', 'off');
                   subplot(1,1,1); megPlotMap(zeros(1,157)); colormap(cmap);
end

% Save sensors of interest that fall within the contour lines
sensorsOfInterest = NaN(size(data));

% Loop over datasets
for ii = 1:size(data,1)
    
    dataToPlot = data(ii,:);
%     colormapLims =  [-1 1]*prctile(dataToPlot, colormapPercentile);
    if ii ==1
        colormapLims = [-40, 40];
    else
        colormapLims = [-2.5, 2.5];
    end

    if ~isempty(contourmapPercentile)
       if contourmapPercentile>10 
           contourLines = [1 1]*prctile(dataToPlot, contourmapPercentile); 
       else
           contourLines = contourmapPercentile; 
       end
    else contourLines = [];
    end

    % Plot data or predictions
    figure(fH1);
    subplot(size(data,1),1,ii);
    [~,ch] = megPlotMap(dataToPlot,colormapLims,fH1,'bipolar',sub_ttl{ii},[],[], ...
        'isolines', contourLines, ...
    ...    'chanindx', dataToPlot > max(contourmapLims), ...
        'pointsymbol', markerType, ... '*'
        'pointsize', 10); hold on;
%     if ii ==1; set(gca, 'CLim', [-50,50]); else set(gca, 'CLim', [-2,2]); end
    % Check if a contour line was requested, if so, save those data and
    % replot in separate figure
    if ~isempty(contourLines)
        c = findobj(gca,'Type','Contour');
%         pp = findobj(gca,'Marker',markerType);
        
        for jj = 1:size(c,1)
            % Change line width
            figure(fH1); 
            c(jj).LineWidth = 4;
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