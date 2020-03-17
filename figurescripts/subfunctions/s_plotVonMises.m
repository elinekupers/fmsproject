% s_plotVonMises

allKappas = pi.*[ 0, logspace(log10(.1),log10(2),10), 100];
thetahat  = 0;
alpha     = linspace(-pi, pi, 1000);


cmap = [linspace(57, 239, 12); linspace(83, 160, 12); linspace(164, 34, 12)];

cmap = cmap'./255;

figure(1); clf; set(gcf, 'Color', 'w'); hold all;


for k = 1:length(allKappas)

    C = 1/(2*pi*besseli(0,allKappas(k)));
    p = C * exp(allKappas(k)*cos(alpha-thetahat));
    
    p  = p ./ sum(p);

 plot(alpha, p, 'Color', cmap(k,:), 'LineWidth', 3);
 
end

set(gca, 'FontSize', 15, 'TickDir', 'out');
xlabel('Phase (rad)')
ylabel('Probability')
xlim([-pi,pi])
legend(cellstr(sprintfc('%1.2f *pi', allKappas/pi)))

legend boxoff
title('Von Mises distributions to sample phase')