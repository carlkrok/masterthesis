function U = ga_MPC(Y0, A_disc,B_disc,prediction_horizon, ...
    numControlVariables, numThrusters, numPropellant, U_lb, U_ub, ...
    Y_lb, Y_ub, Y_lb_dotVec, Y_ub_dotVec, chi_ref, rw_vel_ref)

global mtqData
global rwData
global propulsionData
global simConfig
global satelliteConfiguration


rw_vel = Y0(8:11);

if satelliteConfiguration == 2 && simConfig.enableQuatRef
    attitude_weight = 1e15; % 1e8; % 1e15; % 1e9;
elseif satelliteConfiguration == 1 && simConfig.enableQuatRef
    attitude_weight = 1e15; % 1e8; % 1e15;
else
    attitude_weight = 0;
end
    
attitude_eta_weight = 0;

if simConfig.enableOmegaRef
    omega_weight = 1e15;
else
    omega_weight = 0;
end


if simConfig.enableRW
    rw_actuation_weights = rwData.idlePower + 1/rwData.efficiency .* ...
    rwData.I_mat * (rw_vel_ref .* ones(4,1)); 
    if satelliteConfiguration == 2
        rw_momentum_weight = 10; % 1e-15; % 100
    elseif satelliteConfiguration == 1 
        rw_momentum_weight = 10; % 1e5; % 100
    end
    if simConfig.enableOmegaRef && satelliteConfiguration == 2
        rw_momentum_weight = 100;
    elseif simConfig.enableOmegaRef && satelliteConfiguration == 1 
        rw_momentum_weight = 1e4; % 100
    end
else
    rw_actuation_weights = zeros(4,1);
    rw_momentum_weight = 0;
end

if simConfig.enableMTQ
    mtq_weight = mtqData.powerFactor;
else
    mtq_weight = 0;
end

if simConfig.enablePropulsion
    thrust_weight = propulsionData.power/propulsionData.maxThrust; ...1e3*
    if satelliteConfiguration == 2
        thrust_weight = 1e14*thrust_weight; % 5e8*thrust_weight; % 1e3* % 5e8*thrust_weight
    end
else
    thrust_weight = 0;
end



A_chi = A_disc;

size_b_disc = size(B_disc);
B_chi = [B_disc, zeros(size_b_disc(1),(prediction_horizon-1)*size_b_disc(2))];

U_int_index_0 = 8;
for thrustIter = 2:numThrusters
    U_int_index_0 = [U_int_index_0; 8+thrustIter-1];
end
U_int_index = U_int_index_0;

for horizonIter = 2:prediction_horizon
    
    A_chi = [A_chi; A_disc^horizonIter];
    
    B_chi_row = A_disc^(prediction_horizon-1)*B_disc;
    for bMatIter = 1:(horizonIter-1)
        B_chi_row = [B_chi_row, A_disc^(horizonIter-1-bMatIter)*B_disc];
    end
    B_chi_row = [B_chi_row,zeros(size_b_disc(1),(prediction_horizon-horizonIter)*size_b_disc(2))];
    
    B_chi = [B_chi; B_chi_row];
    
    U_int_index = [U_int_index; U_int_index_0+(horizonIter-1)*numControlVariables];
  
end

chi_0 = A_chi * Y0; % - X_ref;

X_weightMat = diag([attitude_eta_weight,attitude_weight .* ones(1,3), ...
    omega_weight .* ones(1,3), ...
    rw_momentum_weight .* ones(1,4), zeros(1,numPropellant)]);
H_cell = repmat({X_weightMat},1,prediction_horizon);
H_attitude = blkdiag(H_cell{:});

U_weightMat = diag( [mtq_weight .* ones(1,3), ...
    zeros(1,4), thrust_weight .* ones(1,numThrusters)] );
W_cell = repmat({U_weightMat},1,prediction_horizon);
W_actuation = blkdiag(W_cell{:});

Z_actuation = repmat([zeros(1,3), rw_actuation_weights', ...
    zeros(1,numThrusters)],1,prediction_horizon);

A_ub_cell = repmat({diag(Y_ub_dotVec)},1,prediction_horizon);
A_ub = blkdiag(A_ub_cell{:});
A_lb_cell = repmat({-1.*diag(Y_lb_dotVec)},1,prediction_horizon);
A_lb = blkdiag(A_lb_cell{:});

A_y = [ A_ub; A_lb ];

b_y = [repmat(Y_ub,prediction_horizon,1);repmat(-1.*Y_lb,prediction_horizon,1)];

A_u = A_y * B_chi;

b_u = b_y - A_y * chi_0;

c = (( chi_0' * H_attitude * B_chi  - chi_ref'* H_attitude * B_chi ));
D = ( B_chi' * H_attitude * B_chi + W_actuation );
e = Z_actuation;

options = optimoptions('ga', 'Display', 'off', 'MaxTime', 30); %, ...
    %'ConstraintTolerance', 1e-6, 'FunctionTolerance', 1e-9, ...
    %'PopulationSize', 500); %, 'PlotFcn', @gaplotbestf);

x = ga(@( U ) ga_CostFn( U, c, D, e, U_int_index ), ...
    numControlVariables*prediction_horizon,A_u,b_u, ... [],[], ... 
    [],[],repmat(U_lb,prediction_horizon,1), ...
    repmat(U_ub,prediction_horizon,1),[],U_int_index, options);

% options = optimoptions('surrogateopt', 'Display', 'off', 'MaxTime', 5, 'PlotFcn', []);
%     
% x = surrogateopt(@( U ) surrogate_CostFn( U, c, D, e, U_int_index, A_u,b_u ), ...
%     repmat(U_lb,prediction_horizon,1), ...
%     repmat(U_ub,prediction_horizon,1), ...
%     U_int_index, ...
%     options);


%%

U = x(1:numControlVariables)';


end

