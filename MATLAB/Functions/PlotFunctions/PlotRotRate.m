function PlotRotRate( omega, timeVec, omegaRefTime, omegaRef )

figure
hold on
x = timeVec;
tiledlayout(3,1)
ax1 = nexttile;
grid on
plot(x, omega(:,1), omegaRefTime', omegaRef(1,:));%.*(180/pi))
ylabel('[ rad / s ]')
xlabel('[s]')
title('\omega_{sat}^{body} X-Axis')
ax2 = nexttile;
grid on
plot(x, omega(:,2), omegaRefTime', omegaRef(2,:)');%.*(180/pi))
ylabel('[ rad / s ]')
xlabel('[s]')
title('\omega_{sat}^{body} Y-Axis')
ax3 = nexttile;
grid on
plot(x, omega(:,3), omegaRefTime', omegaRef(3,:)');%.*(180/pi))
ylabel('[ rad / s ]')
xlabel('[s]')
title('\omega_{sat}^{body} Z-Axis')
grid(ax1,'on')
grid(ax2,'on')
grid(ax3,'on')
hold off

end

