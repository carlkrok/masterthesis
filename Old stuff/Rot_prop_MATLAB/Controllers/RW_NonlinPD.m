function dw_cmd_rw = RW_NonlinPD( q_error, omega_error, sat_I_mat, ...
    sat_omega, rw_w, rw_A_mat, A_MPinv_mat, rw_I_mat, rw_maxAcc, ...
    rw_maxVel, t_MT, K_p, K_d, RotMat_structToBody )


t_ref_sat_body = - K_p * q_error(2:4) ...
    - K_d * omega_error;
    % - t_MT ...
    
t_ref_sat_struct = RotMat_structToBody' * t_ref_sat_body;   

% vel_ub = rw_maxVel .* ones(length(rw_w), 1) - rw_w;
% vel_lb = -rw_maxVel .* ones(length(rw_w), 1) - rw_w;
% 
% acc_ub = rw_maxAcc .* ones(length(rw_w), 1);
% acc_lb = -rw_maxAcc .* ones(length(rw_w), 1);
% 
% dw_rw_ub = max(min(vel_ub, acc_ub), acc_lb);
% dw_rw_lb = min(max(vel_lb, acc_lb), acc_ub);

t_ref_rw_struct = -t_ref_sat_struct;

% dw_cmd_rw = lsqlin(rw_A_mat * rw_I_mat, t_ref_rw_struct,[],[],[],[], ...
%    dw_rw_lb, ... % Lower bound
%    dw_rw_ub); % Upper bound

dw_cmd_rw = rw_I_mat \ (A_MPinv_mat * t_ref_rw_struct);


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

