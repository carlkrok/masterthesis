function [wdot, t_body] = RW_Torque( w_sat, appliedTorque, ...
    rw_I_mat, h_rw, w_rw_w, A_mat, A_MPinv_mat )

wdot = appliedTorque ./ rw_I_mat(1,1) - A_MPinv_mat * wdot_sat_body;

t_body = -A_mat * appliedTorque;

end
