function PlotDisturbanceTorqueNorm( disturbance_torques_norm, ...
    srp_torque_norm, drag_torque_norm, grav_torque_norm )

figure
hold on
grid on
title('Norm of Disturbance Torques')
plot(disturbance_torques_norm(:,1),disturbance_torques_norm(:,2), '*')
plot(srp_torque_norm(:,1),srp_torque_norm(:,2), '*')
plot(drag_torque_norm(:,1),drag_torque_norm(:,2), '*')
plot(grav_torque_norm(:,1),grav_torque_norm(:,2), '*')
legend('Total','SRP','Drag','Grav')
hold off

end

