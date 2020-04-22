function U = MAP_MPC(Y0, A_disc,B_disc,prediction_horizon, U_lb, U_ub, ...
    Y_lb, Y_ub, Y_lb_dotVec, Y_ub_dotVec, qRef)

global plotData

X_ref = repmat([zeros(7,1);qRef(2:4);zeros(19,1)],prediction_horizon,1);

attitude_weight = 1e8;
rw_actuation_weight = 1;
mtq_weight = 1;
thrust_weight = 1;
rw_momentum_weight = 0.1;


u_integrated_cell = {prediction_horizon.*eye(13)};
for horizonIter = 1:(prediction_horizon-1)
    u_integrated_cell{horizonIter+1} = ...
        (prediction_horizon-horizonIter) .* eye(13);
end


A_chi = A_disc;

size_b_disc = size(B_disc);
B_chi = [B_disc, zeros(size_b_disc(1),(prediction_horizon-1)*size_b_disc(2))];

for horizonIter = 2:prediction_horizon
    
    A_chi = [A_chi; A_disc^horizonIter];
    
    B_chi_row = A_disc^(prediction_horizon-1)*B_disc;
    for bMatIter = 1:(horizonIter-1)
        B_chi_row = [B_chi_row, A_disc^(horizonIter-1-bMatIter)*B_disc];
    end
    B_chi_row = [B_chi_row,zeros(size_b_disc(1),(prediction_horizon-horizonIter)*size_b_disc(2))];
    
    B_chi = [B_chi; B_chi_row];
  
end

chi_0 = A_chi * Y0; % - X_ref;

X_weightMat = diag([zeros(1,7),attitude_weight .* ones(1,3),zeros(1,3), ...
    rw_momentum_weight .* ones(1,4), zeros(1,12)]);
H_cell = repmat({X_weightMat},1,prediction_horizon);
H = blkdiag(H_cell{:});

U_weightMat = diag( [mtq_weight .* ones(1,3), ...
    rw_actuation_weight .* ones(1,4), thrust_weight .* ones(1,6)] );
Y_cell = repmat({U_weightMat},1,prediction_horizon);
Y_weight = blkdiag(Y_cell{:});
U_integratedMat = blkdiag(u_integrated_cell{:});

Y = Y_weight * U_integratedMat;

Ae = eye(length(Y0)*prediction_horizon);
be = chi_0;
pE = B_chi;

A_ub_cell = repmat({diag(Y_ub_dotVec)},1,prediction_horizon);
A_ub = blkdiag(A_ub_cell{:});
A_lb_cell = repmat({diag(Y_lb_dotVec)},1,prediction_horizon);
A_lb = blkdiag(A_lb_cell{:});
A = [ A_ub; A_lb ];

b = [repmat(Y_ub,prediction_horizon,1);repmat(Y_lb,prediction_horizon,1)];

Ath = [eye(length(U_ub)*prediction_horizon); ...
    -1.*eye(length(U_lb)*prediction_horizon)];

bth = [repmat(U_ub,prediction_horizon,1);repmat(-1.*U_lb,prediction_horizon,1)];

opt = Opt( 'H', H, 'Y', Y, 'Ae', Ae, 'be', be, 'pE', pE) %, 'A', A, 'b', b, ...
    % 'verbose', 2);%, 'problem_type', 'MIQP' ) % 'Ath', Ath, 'bth', bth, 

res = opt.solve

%%

  
%[xmin,fmin,flag,Extendedflag]=miqp(H,f,A,b,Aeq,beq,vartype,lb,ub,x0,Options)

cost_attitude = attitude_weight * dot( ...
    ( sol.Y_opt .* Yref_dotVec_tot - Yref_vec_tot ).^2 ...
    , ...
    ones( size( Yref_vec_tot ) ) );
cost_mtq = mtq_weight * dot(u_cont_integrated_mat*sol.u_cont.^2, u_mtq_dotVec);
cost_rw_actuation = rw_actuation_weight * dot(u_cont_integrated_mat*sol.u_cont.^2, u_rw_dotVec);
cost_thrust = thrust_weight * ones(1,6) * u_int_integrated_mat*sol.u_int;
cost_actuation = cost_mtq + cost_rw_actuation + cost_thrust;
cost_rw_momentum = rw_momentum_weight * dot( (sol.Y_opt .* Y_rw_momentum_dotVec_tot).^2, ...
    ones( size( Yref_vec_tot ) ) );


plotData.cost_attitude = [plotData.cost_attitude, cost_attitude];
plotData.cost_actuation = [plotData.cost_actuation, cost_actuation];
plotData.cost_rw_momentum = [plotData.cost_rw_momentum, cost_rw_momentum];

U = [x];


end

