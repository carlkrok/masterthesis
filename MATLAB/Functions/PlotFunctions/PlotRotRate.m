function PlotRotRate( omega )

numSteps = length(omega(:,1));

figure
hold on
x = 0:numSteps-1;
tiledlayout(3,1)
nexttile
grid on
plot(x, omega(:,1).*(180/pi))
title('\omega_{sat}^{body} X-Axis')
nexttile
grid on
plot(x, omega(:,2).*(180/pi))
title('\omega_{sat}^{body} Y-Axis')
nexttile
grid on
plot(x, omega(:,3).*(180/pi))
title('\omega_{sat}^{body} Z-Axis')
hold off

end

