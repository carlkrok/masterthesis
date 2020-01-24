function [ wdot_sat, wdot_rw ] = Forces_RW( t, Y, qRefTraj, sat_I_mat, sat_rw_rotMat )

wdot_rw = zeros(4,1);

torque = zeros(3,1);
wdot_sat = torque;

end

