function PlotMTQDipole( mtq_m )
    numSteps = length(mtq_m(:,1));
    figure
    hold on
    t = mtq_m(:,4);
    tiledlayout(3,1)
    ax1 = nexttile;
    plot(t, mtq_m(:,1))
    title('M^{b} X-Axis')
    ax2 = nexttile;
    plot(t, mtq_m(:,2))
    title('M^{b} Y-Axis')
    ax3 = nexttile;
    plot(t, mtq_m(:,3))
    title('M^{b} Z-Axis')
    grid(ax1,'on')
    grid(ax2,'on')
    grid(ax3,'on')
    hold off
end

