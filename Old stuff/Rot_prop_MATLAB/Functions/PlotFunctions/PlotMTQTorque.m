function PlotMTQTorque( mtq_t )
    numSteps = length(mtq_t(:,1));
    figure
    hold on
    t = mtq_t(:,4);
    tiledlayout(3,1)
    ax1 = nexttile;
    grid on
    plot(t, mtq_t(:,1))
    title('T_{MTQ}^{b} X-Axis')
    ax2 = nexttile;
    grid on
    plot(t, mtq_t(:,2))
    title('T_{MTQ}^{b} Y-Axis')
    ax3 = nexttile;
    grid on
    plot(t, mtq_t(:,3))
    title('T_{MTQ}^{b} Z-Axis')
    grid(ax1,'on')
    grid(ax2,'on')
    grid(ax3,'on')
    hold off
end

