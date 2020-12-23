function [] = visualizePosteriorSensors1D(dataToPlot, fig_ttl, sub_ttl, saveFigures, figureDir)

% get sensor locations
load(which('meg160_example_hdr.mat'))
layout = ft_prepare_layout([],hdr);
xpos   = layout.pos(1:157,1);
ypos   = layout.pos(1:157,2);

xbins = linspace(min(xpos),max(xpos),11);
for ii = 1:length(xbins)-1
        xposbin = ((xpos>=xbins(ii)) & (xpos<xbins(ii+1)));
        posteriorSensorLoc(ii,:) = (ypos<0 & xposbin);
        
        mn1Ddata(1,ii) = nanmean(dataToPlot(1,posteriorSensorLoc(ii,:)),2);
        mn1Ddata(2,ii) = nanmean(dataToPlot(2,posteriorSensorLoc(ii,:)),2);
        
end

xq = linspace(1,10,20);
splineData(1,:) = spline(1:10,mn1Ddata(1,:),xq);
splineData(2,:) = spline(1:10,mn1Ddata(2,:),xq);

x = (xq-1);
x = x-(max(x)/2);

figure; set(gcf, 'color', 'w');
subplot(211);
plot(x,splineData(1,:), 'r-','LineWidth',2); hold on;
if any(splineData(1,:)<0)
    plot(x,zeros(size(x)), 'k','LineWidth',1);
end
box off; title(sub_ttl{1})
set(gca, 'TickDir','out','ticklength',[0.010 0.010], 'FontSize',15);
ylabel('Amplitude (fT)')
subplot(212);
plot(x,splineData(2,:), 'b-','LineWidth',2); hold on;
if any(splineData(2,:)<0)
    plot(x,zeros(size(x)), 'k','LineWidth',1);
end
box off; title(sub_ttl{1})
set(gca, 'TickDir','out','ticklength',[0.010 0.010],'FontSize',15);
xlabel('x position sensor map'); ylabel('Power (fT^2)')

% Save if requested
if saveFigures
    figurewrite(fullfile(figureDir,fig_ttl),[],0,'.',1);
end

