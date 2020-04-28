function PlotCosts(cost_attitude,cost_actuation,cost_rw_momentum)

figure
hold on
grid on
title('MPC Costs')
time = 1:length(cost_attitude);
plot(time,cost_attitude)
plot(time,cost_actuation)
plot(time,cost_rw_momentum)
legend('Attitude','Actuation','RW Momentum')
hold off

end

