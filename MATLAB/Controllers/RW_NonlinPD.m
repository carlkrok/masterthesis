function t_ref_RW = RW_NonlinPD( q_error, omega_error, sat_I_mat, ...
    sat_omega, rw_momentum, rw_A_mat, rw_AInv_mat, rw_maxTorque, ...
    rw_maxMomentum, t_MT, controller_rw_natFreq, ...
    controller_rw_dampRatio )

t_ref_sat = -2 * controller_rw_natFreq^2 * sat_I_mat * q_error(2:4) ...
    - 2 * controller_rw_dampRatio * controller_rw_natFreq * sat_I_mat ...
    * omega_error ...
    - t_MT ...
    + cross( sat_omega, (sat_I_mat * sat_omega + rw_A_mat * rw_momentum )); 

momentum_ub = rw_maxMomentum .* ones(length(rw_momentum), 1) - rw_momentum;
momentum_lb = -rw_maxMomentum .* ones(length(rw_momentum), 1) - rw_momentum;

torque_ub = rw_maxTorque .* ones(length(rw_momentum), 1);
torque_lb = -rw_maxTorque .* ones(length(rw_momentum), 1);

t_rw_ub = max(min(momentum_ub, torque_ub), torque_lb);
t_rw_lb = min(max(momentum_lb, torque_lb), torque_ub);

if any(t_rw_lb > t_rw_ub, 'all')
    disp("Optimization problem")
end

t_ref_RW = lsqlin(rw_A_mat,t_ref_sat,[],[],[],[], ...
    t_rw_lb, ... % Lower bound
    t_rw_ub); % Upper bound

if any(isnan(t_ref_RW), 'all')
    error("NaN found in RW controller")
end

end

