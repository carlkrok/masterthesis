function [hdot, t_body] = RW_Torque( w_sat, appliedTorque, h_rw, A_mat, A_MPinv_mat )

hdot = appliedTorque;

t_body = -A_mat * appliedTorque - ( cross( w_sat, A_mat * h_rw ) );

end
