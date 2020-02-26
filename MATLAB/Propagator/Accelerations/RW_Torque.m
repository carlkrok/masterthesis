function [hdot, t_body] = RW_Torque( w_sat, appliedTorque, h_rw, A_mat, A_MPinv_mat )

hdot = A_MPinv_mat * ( cross( A_mat * h_rw, w_sat ) ) - appliedTorque;

t_body = -1 * A_mat * hdot - cross( w_sat, A_mat * h_rw );

end
