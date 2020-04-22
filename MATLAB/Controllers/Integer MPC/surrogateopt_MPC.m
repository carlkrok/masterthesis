function U = surrogateopt_MPC(Y0, A_disc,B_disc,prediction_horizon, U_lb, U_ub, ...
    Y_lb, Y_ub, Y_lb_dotVec, Y_ub_dotVec, qRef)

global plotData

attitude_weight = 1e5;
rw_actuation_weight = 10;
mtq_weight = 10;
thrust_weight = 1e-10;
rw_momentum_weight = 0.1;

X_ref = repmat([zeros(7,1);qRef(2:4);zeros(19,1)],prediction_horizon,1);

U_lb_tot = repmat(U_lb,prediction_horizon,1);
U_ub_tot = repmat(U_ub,prediction_horizon,1);

Y_lb_tot = repmat(Y_lb,prediction_horizon,1);
Y_ub_tot = repmat(Y_ub,prediction_horizon,1);

Y_lb_dotVec_tot = repmat(Y_lb_dotVec,prediction_horizon,1);
Y_ub_dotVec_tot = repmat(Y_ub_dotVec,prediction_horizon,1);

A_chi = A_disc;

B_disc = B_disc;
size_b_disc = size(B_disc);
B_chi = [B_disc, zeros(size_b_disc(1),(prediction_horizon-1)*size_b_disc(2))];

U_int_index_0 = [8, 9, 10, 11, 12, 13]';
U_int_index = U_int_index_0;

for horizonIter = 2:prediction_horizon
    
    A_chi = [A_chi; A_disc^horizonIter];
    
    B_chi_row = A_disc^(prediction_horizon-1)*B_disc;
    for bMatIter = 1:(horizonIter-1)
        B_chi_row = [B_chi_row, A_disc^(horizonIter-1-bMatIter)*B_disc];
    end
    B_chi_row = [B_chi_row,zeros(size_b_disc(1),(prediction_horizon-horizonIter)*size_b_disc(2))];
    
    B_chi = [B_chi; B_chi_row];
    
    U_int_index = [U_int_index; U_int_index_0+(horizonIter-1)*13];
  
end

chi_0 = A_chi * Y0; % - X_ref;

u_integrated_cell = {prediction_horizon.*eye(13)};
for horizonIter = 1:(prediction_horizon-1)
    u_integrated_cell{horizonIter+1} = ...
        (prediction_horizon-horizonIter) .* eye(13);
end

X_weightMat = diag([zeros(1,7),attitude_weight .* ones(1,3),zeros(1,3), ...
    rw_momentum_weight .* ones(1,4), zeros(1,12)]);
H_cell = repmat({X_weightMat},1,prediction_horizon);
H = blkdiag(H_cell{:});

U_weightMat = diag( [mtq_weight .* ones(1,3), ...
    rw_actuation_weight .* ones(1,4), thrust_weight .* ones(1,6)] );
W_cell = repmat({U_weightMat},1,prediction_horizon);
W_weight = blkdiag(W_cell{:});
U_integratedMat = blkdiag(u_integrated_cell{:});

W = W_weight * U_integratedMat;
  
options = optimoptions(@surrogateopt, 'Display', 'off', ...
    'MaxFunctionEvaluations', 200, 'UseParallel', false, ...
    'MinSampleDistance', 1e-3);
[sol, fval] = surrogateopt(@(U) MPC_costfn( U, Y_lb_tot, Y_ub_tot, ...
    Y_lb_dotVec_tot, Y_ub_dotVec_tot, B_chi, chi_0, H, W ), ...
    U_lb_tot, U_ub_tot, U_int_index, options);

U = sol(1:13)';


end

