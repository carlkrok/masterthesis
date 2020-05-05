function PlotMTQDipole( mtq_m )
figure
hold on
t = mtq_m(:,1);
tiledlayout(3,1)
ax1 = nexttile;
plot(t, mtq_m(:,2))
title('M^{b} X-Axis')
ylabel('[Am^2]')
xlabel('[s]')
ax2 = nexttile;
plot(t, mtq_m(:,3))
title('M^{b} Y-Axis')
ylabel('[Am^2]')
xlabel('[s]')
ax3 = nexttile;
plot(t, mtq_m(:,4))
title('M^{b} Z-Axis')
ylabel('[Am^2]')
xlabel('[s]')
grid(ax1,'on')
grid(ax2,'on')
grid(ax3,'on')
hold off
end

