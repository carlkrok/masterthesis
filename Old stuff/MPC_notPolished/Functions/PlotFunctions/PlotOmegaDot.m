function PlotOmegaDot( wdot_sat_body )

figure
hold on
x = wdot_sat_body(:,1);
tiledlayout(3,1)
ax1 = nexttile;
grid on
plot(x, wdot_sat_body(:,2),'*');%.*(180/pi))
ylabel('rad/s')
title('d^b/dt \omega_{sat}^{body} X-Axis')
ax2 = nexttile;
grid on
plot(x, wdot_sat_body(:,3),'*');%.*(180/pi))
ylabel('rad/s')
title('d^b/dt \omega_{sat}^{body} Y-Axis')
ax3 = nexttile;
grid on
plot(x, wdot_sat_body(:,4),'*');%.*(180/pi))
ylabel('rad/s')
title('d^b/dt \omega_{sat}^{body} Z-Axis')
grid(ax1,'on')
grid(ax2,'on')
grid(ax3,'on')
hold off

end
