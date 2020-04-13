function t_cmd_rw = RW_NonlinPD( q_error, omega_error, sat_I_mat, ...
    sat_omega, rw_w, rw_A_mat, rw_I_mat, rw_maxAcc, ...
    rw_maxVel, t_MT, K_p, K_d )

% K_p = 1;
% K_d = 0.1;
% t_ref_sat =  - K_p * q_error(2:4) + K_d * omega_error;

rw_maxTorque = rw_I_mat(1,1) * rw_maxAcc;
rw_maxMomentum = rw_I_mat(1,1) * rw_maxVel;

t_ref_sat = - K_p * sat_I_mat * q_error(2:4) ...
    - K_d * sat_I_mat * omega_error ...
    + cross( sat_omega, (sat_I_mat * sat_omega + rw_A_mat * rw_w )); 
    % - t_MT ...

momentum_ub = rw_maxMomentum .* ones(length(rw_w), 1) - rw_w;
momentum_lb = -rw_maxMomentum .* ones(length(rw_w), 1) - rw_w;

torque_ub = rw_maxTorque .* ones(length(rw_w), 1);
torque_lb = -rw_maxTorque .* ones(length(rw_w), 1);

t_rw_ub = max(min(momentum_ub, torque_ub), torque_lb);
t_rw_lb = min(max(momentum_lb, torque_lb), torque_ub);

t_ref_rw = -t_ref_sat;

t_cmd_rw = lsqlin(rw_A_mat,t_ref_rw,[],[],[],[], ...
    t_rw_lb, ... % Lower bound
    t_rw_ub); % Upper bound

if any(isnan(t_cmd_rw), 'all')
    pause("NaN found in RW controller")
end

end
