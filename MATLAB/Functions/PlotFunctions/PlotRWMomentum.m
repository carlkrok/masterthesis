function PlotRWMomentum( rw_momentum )

numSteps = length(rw_momentum(:,1));

figure
hold on
grid on
x = 0:numSteps-1;
tiledlayout(4,1)
ax1 = nexttile;
plot(x, rw_momentum(:,1) .* (30 / pi))
ylabel('RPM')
title('RW 1 velocity')
ax2 = nexttile;
plot(x, rw_momentum(:,2) .* (30 / pi))
ylabel('RPM')
title('RW 2 velocity')
ax3 = nexttile;
plot(x, rw_momentum(:,3) .* (30 / pi))
ylabel('RPM')
title('RW 3 velocity')
ax4 = nexttile;
plot(x, rw_momentum(:,4) .* (30 / pi))
ylabel('RPM')
title('RW 4 velocity')
grid(ax1,'on')
grid(ax2,'on')
grid(ax3,'on')
grid(ax4,'on')
hold off

end

