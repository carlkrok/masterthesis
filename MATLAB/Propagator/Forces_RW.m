function [ moment_body ] = Forces_RW( w_rw, rw_A_mat, rw_I_mat, rw_maxTorque, rw_maxMomentum )

moment_body = rw_A_mat*rw_I_mat*w_rw;

end

