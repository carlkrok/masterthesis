function [hdot, t_body] = RW_Torque( w_sat, appliedTorque, h_rw, A_mat, A_MPinv_mat )

hdot = -1 .* A_MPinv_mat * ( cross( w_sat, A_mat * h_rw ) ) + appliedTorque;

t_rw = A_mat * hdot + cross( w_sat, A_mat * h_rw );

t_body = -t_rw;

end
