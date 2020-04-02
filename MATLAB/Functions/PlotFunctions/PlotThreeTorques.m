function PlotThreeTorques( t1, t2, t3 )
    figure
    hold on
    time1 = t1(:,4);
    time2 = t2(:,4);
    time3 = t3(:,4);
    tiledlayout(3,1)
    ax1 = nexttile;
    grid on
    plot(time1, t1(:,1),'r',time2, t2(:,1),'g',time3, t3(:,1),'b')
    title('Torques X-Axis')
    ax2 = nexttile;
    grid on
    plot(time1, t1(:,2),'r',time2, t2(:,2),'g',time3, t3(:,2),'b')
    title('Torques Y-Axis')
    ax3 = nexttile;
    grid on
    plot(time1, t1(:,3),'r',time2, t2(:,3),'g',time3, t3(:,3),'b')
    title('Torques Z-Axis')
    grid(ax1,'on')
    grid(ax2,'on')
    grid(ax3,'on')
    hold off
end