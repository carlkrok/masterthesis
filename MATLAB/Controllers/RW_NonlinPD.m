function dw_cmd_rw = RW_NonlinPD( q_error, omega_error, sat_I_mat, ...
    rw_w, A_MPinv_mat, rw_I_mat, rw_maxAcc, ...
    rw_maxVel, t_MT, K_p, K_d, K_dd, acc_error )


t_ref_sat_body = - K_p * sat_I_mat * q_error(2:4) ...
    - K_d * sat_I_mat * omega_error ...
    - K_dd * sat_I_mat * acc_error;
    % - t_MT ... 
      

t_ref_rw = -t_ref_sat_body;

dw_cmd_rw = rw_I_mat \ (A_MPinv_mat * t_ref_rw);


for rwIter = 1:length(dw_cmd_rw)
    if dw_cmd_rw(rwIter) > rw_maxAcc
        dw_cmd_rw(rwIter) = rw_maxAcc;
    elseif dw_cmd_rw(rwIter) < -rw_maxAcc
            dw_cmd_rw(rwIter) = -rw_maxAcc;
    end
    
    if rw_w(rwIter) >= rw_maxVel && dw_cmd_rw(rwIter) > 0
        rw_w(rwIter) = rw_maxVel;
        dw_cmd_rw(rwIter) = -1e-3;
    elseif rw_w(rwIter) <= -rw_maxVel && dw_cmd_rw(rwIter) < 0
        rw_w(rwIter) = -rw_maxVel;
        dw_cmd_rw(rwIter) = 1e-3;
    end
end


end

