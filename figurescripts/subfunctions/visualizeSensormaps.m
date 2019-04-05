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

if ~isempty(contourmapPercentile) && length(contourmapPercentile)==1
    fH2 = figure; clf; set(fH2,'position',[ 1400, 923, 700, 473], 'Name', fig_ttl{2}, 'NumberTitle', 'off');
                   subplot(1,1,1); megPlotMap(zeros(1,157)); colormap(cmap);
end

% Save sensors of interest that fall within the contour lines
sensorsOfInterest = NaN(size(data));

% Loop over datasets
for ii = 1:size(data,1)
    
    dataToPlot = data(ii,:);
    colormapLims =  [-1 1]*prctile(dataToPlot, colormapPercentile);
    
    if ~isempty(contourmapPercentile) && length(contourmapPercentile)==1      
        contourmapLims = [1 1]*prctile(dataToPlot, contourmapPercentile);
    else contourmapLims = [];
    end

    % Plot data or predictions
    figure(fH1);
    subplot(size(data,1),1,ii);
    [~,ch] = megPlotMap(dataToPlot,colormapLims,fH1,'bipolar',[],[],[], ...
        'isolines', contourmapLims, ...
    ...    'chanindx', dataToPlot > max(contourmapLims), ...
        'pointsymbol', markerType, ... '*'
        'pointsize', 10); hold on;
    
    % Check if a contour line was requested, if so, save those data
    if ~isempty(contourmapPercentile) && length(contourmapPercentile)==1
        c = findobj(gca,'Type','Contour'); c.LineWidth = 4;
        pp = findobj(gca,'Marker',markerType);
        
    % Check if a contour lines are a matrix, if so, plot contour lines based on these data   
    elseif ~isempty(contourmapPercentile) && length(contourmapPercentile)>1        
        
        fHP = figure(99); clf; subplot(3,1,1);
        external_contourmapLims1 = [1 1]*prctile(contourmapPercentile(idx,:), 93.6);
        external_colormapLims1   = [-1 1]*prctile(contourmapPercentile(idx,:), colormapPercentile);
        megPlotMap(contourmapPercentile(idx,:),external_colormapLims1,[],bipolar,[],[],[],'isolines', external_contourmapLims1);
        cP1 = findobj(gca,'Type','Contour');
        
        subplot(3,1,2);
        external_contourmapLims2 = [1 1]*prctile(contourmapPercentile(idx+1,:), 93.6);
        external_colormapLims2   = [-1 1]*prctile(contourmapPercentile(idx+1,:), colormapPercentile);
        megPlotMap(contourmapPercentile(idx+1,:),external_colormapLims2,[],bipolar,[],[],[],'isolines', external_contourmapLims2);
        cP2 = findobj(gca,'Type','Contour');
        
        subplot(3,1,3);
        external_contourmapLims3 = [1 1]*prctile(contourmapPercentile(idx+2,:), 93.6);
        external_colormapLims3   = [-1 1]*prctile(contourmapPercentile(idx+2,:), colormapPercentile);
        megPlotMap(contourmapPercentile(idx+2,:),external_colormapLims3,[],bipolar,[],[],[],'isolines', external_contourmapLims3);
        cP3 = findobj(gca,'Type','Contour');
        
         % Plot them in actual mesh
        set(0, 'currentfigure', fH1);
        contour(cP1.XData, cP1.YData, cP1.ZData, external_contourmapLims1, 'LineColor','k', 'Fill','off','LineWidth',2);
        contour(cP2.XData, cP2.YData, cP2.ZData, external_contourmapLims2, 'LineColor','w', 'Fill','off','LineWidth',2);
        contour(cP3.XData, cP3.YData, cP3.ZData, external_contourmapLims3, 'LineColor','w', 'Fill','off','LineWidth',2);

        sensorsOfInterest(idx,:) = contourmapPercentile(idx,:) > max(external_contourmapLims1);
        sensorsOfInterest(idx+1,:) = contourmapPercentile(idx+1,:) > max(external_contourmapLims2);     
        sensorsOfInterest(idx+2,:) = contourmapPercentile(idx+2,:) > max(external_contourmapLims3);        

    end
    
    % Make figure pretty
    set(ch,'box','off','tickdir','out','ticklength',[0.010 0.010], 'FontSize',12); title(sub_ttl{ii})
        
    % Plot overlap
    if ~isempty(contourmapLims)
        figure(fH2); hold all;   
        contourf(c.XData, c.YData, c.ZData, contourmapLims, 'LineColor',colorMarkers{ii}, 'Fill','off','LineWidth',4);
        %%scatter(pp(1).XData,pp(1).YData, 150, colorMarkers{ii},'*'); 
        colorbar off;
    
        sensorsOfInterest(ii,:) = dataToPlot > max(contourmapLims);
    end
    
end

% Save if requested
if saveFigures
    if ~exist('figureDir', 'dir'); mkdir(figureDir); end
    set(0, 'currentfigure', fH1);
    figurewrite(fullfile(figureDir,fig_ttl{1}),[],0,'.',1);
    if ~isempty(contourmapLims)
        set(0, 'currentfigure', fH2);
        figurewrite(fullfile(figureDir,fig_ttl{2}),[],0,'.',1);
    end
end