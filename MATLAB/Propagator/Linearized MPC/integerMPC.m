function U = integerMPC(Y0, A_disc,B_disc,prediction_horizon, U_lb, U_ub, ...
    Y_lb, Y_ub, Y_lb_dotVec, Y_ub_dotVec, qRef)

Yref_vec = [zeros(7,1);qRef(2:4);zeros(19,1)];
Yref_vec_tot = repmat(Yref_vec,prediction_horizon,1);

Yref_dotVec = [zeros(7,1);ones(3,1);zeros(19,1)];
Yref_dotVec_tot = repmat(Yref_dotVec,prediction_horizon,1);

u_mtq_dotVec = [ones(3,1);zeros(4,1)];
u_rw_dotVec = [zeros(3,1);ones(4,1)];
u_thrust_dotVec = ones(6,1);

U_lb_cont_tot = repmat(U_lb(1:7),prediction_horizon,1);
U_ub_cont_tot = repmat(U_ub(1:7),prediction_horizon,1);
U_lb_int_tot = repmat(U_lb(8:13),prediction_horizon,1);
U_ub_int_tot = repmat(U_ub(8:13),prediction_horizon,1);

Y_lb_tot = repmat(Y_lb,prediction_horizon,1);
Y_ub_tot = repmat(Y_ub,prediction_horizon,1);

Y_lb_dotVec_tot = repmat(Y_lb_dotVec,prediction_horizon,1);
Y_ub_dotVec_tot = repmat(Y_ub_dotVec,prediction_horizon,1);

A_chi = A_disc;

B_disc_cont = B_disc(:,1:7);
size_b_disc_cont = size(B_disc_cont);
B_chi_cont = [B_disc_cont, zeros(size_b_disc_cont(1),(prediction_horizon-1)*size_b_disc_cont(2))];

B_disc_int = B_disc(:,8:13);
size_b_disc_int = size(B_disc_int);
B_chi_int = [B_disc_int, zeros(size_b_disc_int(1),(prediction_horizon-1)*size_b_disc_int(2))];

for horizonIter = 2:prediction_horizon
    
    A_chi = [A_chi; A_disc^horizonIter];
    
    B_chi_cont_row = A_disc^(prediction_horizon-1)*B_disc_cont;
    B_chi_int_row = A_disc^(prediction_horizon-1)*B_disc_int;
    for bMatIter = 1:(horizonIter-1)
        B_chi_cont_row = [B_chi_cont_row, A_disc^(horizonIter-1-bMatIter)*B_disc_cont];
        B_chi_int_row = [B_chi_int_row, A_disc^(horizonIter-1-bMatIter)*B_disc_int];
    end
    B_chi_cont_row = [B_chi_cont_row,zeros(size_b_disc_cont(1),(prediction_horizon-horizonIter)*size_b_disc_cont(2))];
    B_chi_int_row = [B_chi_int_row,zeros(size_b_disc_int(1),(prediction_horizon-horizonIter)*size_b_disc_int(2))];
    
    B_chi_cont = [B_chi_cont; B_chi_cont_row];
    B_chi_int = [B_chi_int; B_chi_int_row];
           
end

chi_0 = A_chi * Y0;

prob = optimproblem('ObjectiveSense','min');

u_delta_int = optimvar('u_int',6*prediction_horizon,1,'Type','integer',...
    'LowerBound',repmat(zeros(6,1),prediction_horizon,1),'UpperBound',U_ub_int_tot);
u_delta_cont = optimvar('u_cont',7*prediction_horizon,1,...
    'LowerBound',U_lb_cont_tot,'UpperBound',U_ub_cont_tot);
Y = optimvar('Y',length(Yref_dotVec_tot),1);

u_cont_integrated_mat = prediction_horizon.*eye(7);
u_int_integrated_mat = prediction_horizon.*eye(6);
for horizonIter = 1:(prediction_horizon-1)
    u_cont_integrated_mat = [ ...
        u_cont_integrated_mat, (prediction_horizon-horizonIter) .* eye(7)];
    u_int_integrated_mat = [ ...
        u_int_integrated_mat, (prediction_horizon-horizonIter) .* eye(6)];
end

attitude_weight = 10;
rw_weight = 1;
mtq_weight = 1;
thrust_weight = 1;

prob.Objective = attitude_weight * dot( ...
    ( Y .* Yref_dotVec_tot - Yref_vec_tot ).^2 ...
    , ...
    ones( size( Yref_vec_tot ) ) ...
    ) ...
    + mtq_weight * dot((u_cont_integrated_mat*u_delta_cont).^2, u_mtq_dotVec) ...
    + rw_weight * dot((u_cont_integrated_mat*u_delta_cont).^2, u_rw_dotVec) ...
    + thrust_weight * dot((u_int_integrated_mat*u_delta_int).^2, u_thrust_dotVec);

cons1 = Y == chi_0 + B_chi_int * u_delta_int + B_chi_cont * u_delta_cont;
prob.Constraints.cons1 = cons1;

cons2 = Y .* Y_lb_dotVec_tot >= Y_lb_tot;
prob.Constraints.cons4 = cons2;

cons3 = Y .* Y_ub_dotVec_tot <= Y_ub_tot;
prob.Constraints.cons5 = cons3;

%show(prob);

options = optimoptions('Display','off');
sol = solve(prob,'options',options);

cost_attitude = attitude_weight * dot( ...
    ( sol.Y .* Yref_dotVec_tot - Yref_vec_tot ).^2 ...
    , ...
    ones( size( Yref_vec_tot ) ) )

cost_actuation = mtq_weight * dot(u_cont_integrated_mat*sol.u.^2, u_mtq_dotVec) ...
    + rw_weight * dot(u_cont_integrated_mat*sol.u.^2, u_rw_dotVec) ...
    + thrust_weight * dot(u_cont_integrated_mat*sol.u.^2, u_thrust_dotVec)

cost = cost_attitude + cost_actuation;

U = sol.u(1:13);


end

