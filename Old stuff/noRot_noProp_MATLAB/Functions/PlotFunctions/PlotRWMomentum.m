function PlotRWMomentum( rw_momentum )

numSteps = length(rw_momentum(:,1));

figure
hold on
grid on
x = 0:numSteps-1;
tiledlayout(4,1)
ax1 = nexttile;
plot(x, rw_momentum(:,1))
title('RW 1 momentum')
ax2 = nexttile;
plot(x, rw_momentum(:,2))
title('RW 2 momentum')
ax3 = nexttile;
plot(x, rw_momentum(:,3))
title('RW 3 momentum')
ax4 = nexttile;
plot(x, rw_momentum(:,4))
title('RW 4 momentum')
grid(ax1,'on')
grid(ax2,'on')
grid(ax3,'on')
grid(ax4,'on')
hold off

end

