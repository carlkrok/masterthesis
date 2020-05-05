function PlotRWMomentum( rw_momentum, timeVec )

figure
hold on
grid on
x = timeVec;
tiledlayout(4,1)
ax1 = nexttile;
plot(x, rw_momentum(:,1) .* (30 / pi))
ylabel('[RPM]')
xlabel('[s]')
title('RW 1 Velocity')
ax2 = nexttile;
plot(x, rw_momentum(:,2) .* (30 / pi))
ylabel('[RPM]')
xlabel('[s]')
title('RW 2 Velocity')
ax3 = nexttile;
plot(x, rw_momentum(:,3) .* (30 / pi))
ylabel('[RPM]')
xlabel('[s]')
title('RW 3 Velocity')
ax4 = nexttile;
plot(x, rw_momentum(:,4) .* (30 / pi))
ylabel('[RPM]')
xlabel('[s]')
title('RW 4 Velocity')
grid(ax1,'on')
grid(ax2,'on')
grid(ax3,'on')
grid(ax4,'on')
hold off

end

