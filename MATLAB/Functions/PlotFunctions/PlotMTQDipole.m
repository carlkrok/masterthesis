function PlotMTQDipole( mtq_m )
figure
hold on
t = extendVecToIntegerTime(mtq_m(:,1));
t = t(2:end);
tiledlayout(3,1)
ax1 = nexttile;
mtq1 = extendVecToIntegerTime(mtq_m(:,2));
mtq1 = mtq1(1:end-1);
plot(t, mtq1)
title('M^{b} X-Axis')
ylabel('[Am^2]')
xlabel('[s]')
ax2 = nexttile;
mtq2 = extendVecToIntegerTime(mtq_m(:,3));
mtq2 = mtq2(1:end-1);
plot(t, mtq2)
title('M^{b} Y-Axis')
ylabel('[Am^2]')
xlabel('[s]')
ax3 = nexttile;
mtq3 = extendVecToIntegerTime(mtq_m(:,4));
mtq3 = mtq3(1:end-1);
plot(t, mtq3)
title('M^{b} Z-Axis')
ylabel('[Am^2]')
xlabel('[s]')
grid(ax1,'on')
grid(ax2,'on')
grid(ax3,'on')
hold off
end

