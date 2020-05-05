function PlotCoM( com, timeVec )

figure
hold on
time = timeVec;
tiledlayout(3,1)
ax1 = nexttile;
grid on
plot(time, com(:,1))
title('CoM X-Axis')
ax2 = nexttile;
grid on
plot(time, com(:,2))
title('CoM Y-Axis')
ax3 = nexttile;
grid on
plot(time, com(:,3))
title('CoM Z-Axis')
grid(ax1,'on')
grid(ax2,'on')
grid(ax3,'on')
hold off

end

