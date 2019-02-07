% s_plotVonMises

allKappas = sort([linspace(0, 5*pi, 21), 0.125*pi, 10*pi]);
thetahat = 0;
alpha = linspace(-pi, pi, 1000);

cmap = parula(length(allKappas));
figure(1); clf; set(gcf, 'Color', 'w'); hold all;


for k = 1:length(allKappas)

    C = 1/(2*pi*besseli(0,allKappas(k)));
    p = C * exp(allKappas(k)*cos(alpha-thetahat));


 plot(alpha, p, 'Color', cmap(k,:), 'LineWidth', 3);
 
end

set(gca, 'FontSize', 15, 'TickDir', 'out');
xlabel('Phase (rad)')
ylabel('Probability')
xlim([-pi,pi])
legend(cellstr(sprintfc('%1.2f *pi', allKappas/pi)))

legend boxoff
title('Von Mises distributions to sample phase')