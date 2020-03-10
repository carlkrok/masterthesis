function PlotRotRate( omega )

numSteps = length(omega(:,1));

figure
hold on
x = 0:numSteps-1;
tiledlayout(3,1)
ax1 = nexttile;
grid on
plot(x, omega(:,1).*(180/pi))
title('\omega_{sat}^{body} X-Axis')
ax2 = nexttile;
grid on
plot(x, omega(:,2).*(180/pi))
title('\omega_{sat}^{body} Y-Axis')
ax3 = nexttile;
grid on
plot(x, omega(:,3).*(180/pi))
title('\omega_{sat}^{body} Z-Axis')
grid(ax1,'on')
grid(ax2,'on')
grid(ax3,'on')
hold off

end

