function U = customMPC(Y0, A_disc,B_disc,prediction_horizon, U_lb, U_ub, ...
    Y_lb, Y_ub, Y_lb_dotVec, Y_ub_dotVec, qRef)

Yref_vec = [zeros(7,1);qRef(2:4);zeros(19,1)];
Yref_vec_tot = repmat(Yref_vec,prediction_horizon,1);

Yref_dotVec = [zeros(7,1);ones(3,1);zeros(19,1)];
Yref_dotVec_tot = repmat(Yref_dotVec,prediction_horizon,1);

u_mtq_dotVec = [ones(3,1);zeros(10,1)];
u_rw_dotVec = [zeros(3,1);ones(4,1);zeros(6,1)];
u_thrust_dotVec = [zeros(7,1);ones(6,1)];
% u_mtq_dotVec_tot = repmat(u_mtq_dotVec,prediction_horizon,1);
% u_rw_dotVec_tot = repmat(u_rw_dotVec,prediction_horizon,1);
% u_thrust_dotVec_tot = repmat(u_thrust_dotVec,prediction_horizon,1);

U_lb_tot = repmat(U_lb,prediction_horizon,1);
U_ub_tot = repmat(U_ub,prediction_horizon,1);

Y_lb_tot = repmat(Y_lb,prediction_horizon,1);
Y_ub_tot = repmat(Y_ub,prediction_horizon,1);

Y_lb_dotVec_tot = repmat(Y_lb_dotVec,prediction_horizon,1);
Y_ub_dotVec_tot = repmat(Y_ub_dotVec,prediction_horizon,1);

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

chi_0 = A_chi * Y0;

prob = optimproblem('ObjectiveSense','min');

%x = optimvar('x',3,3,3,'Type','integer','LowerBound',0,'UpperBound',1)
u_delta = optimvar('u',length(U_lb)*prediction_horizon); %,1,'LowerBound',U_lb_tot,'UpperBound',U_ub_tot);
Y = optimvar('Y',length(Yref_dotVec_tot),1);

u_integrated_mat = [ ...
    repmat(eye(length(U_lb)),1,prediction_horizon); ...
    ];

attitude_weight = 10;
rw_weight = 1;
mtq_weight = 1;
thrust_weight = 1;

prob.Objective = attitude_weight * dot( ...
    ( Y .* Yref_dotVec_tot - Yref_vec_tot ).^2 ...
    , ...
    ones( size( Yref_vec_tot ) ) ...
    ) ...
    + mtq_weight * dot(u_integrated_mat*u_delta.^2, u_mtq_dotVec) ...
    + rw_weight * dot(u_integrated_mat*u_delta.^2, u_rw_dotVec) ...
    + thrust_weight * dot(u_integrated_mat*u_delta.^2, u_thrust_dotVec);

cons1 = Y == chi_0 + B_chi * u_delta;
prob.Constraints.cons1 = cons1;

cons2 = u_delta <= U_ub_tot;
prob.Constraints.cons2 = cons2;

cons3 = u_delta >= U_lb_tot;
prob.Constraints.cons3 = cons3;

cons4 = Y .* Y_lb_dotVec_tot >= Y_lb_tot;
prob.Constraints.cons4 = cons4;

cons5 = Y .* Y_ub_dotVec_tot <= Y_ub_tot;
prob.Constraints.cons5 = cons5;

%show(prob);

options = optimoptions(@quadprog,'Display','off');
sol = solve(prob,'options',options);

cost_attitude = attitude_weight * dot( ...
    ( sol.Y .* Yref_dotVec_tot - Yref_vec_tot ).^2 ...
    , ...
    ones( size( Yref_vec_tot ) ) )

cost_actuation = mtq_weight * dot(u_integrated_mat*sol.u.^2, u_mtq_dotVec) ...
    + rw_weight * dot(u_integrated_mat*sol.u.^2, u_rw_dotVec) ...
    + thrust_weight * dot(u_integrated_mat*sol.u.^2, u_thrust_dotVec)

cost = cost_attitude + cost_actuation;

U = sol.u(1:13);


end

