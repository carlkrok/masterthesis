function U = MAP_MPC(Y0, A_disc,B_disc,prediction_horizon, U_lb, U_ub, ...
    Y_lb, Y_ub, Y_lb_dotVec, Y_ub_dotVec, qRef)

global plotData


temp_Y = repmat(Y0, prediction_horizon, 1);
temp_u_int = repmat(zeros(6,1), prediction_horizon, 1);
temp_u_cont = repmat(zeros(7,1), prediction_horizon, 1);

delta_u_cont = repmat([0.1;0.1;0.1;1;1;1;1],prediction_horizon,1);

cost_total = 0;

attitude_weight = 1e6;
rw_actuation_weight = 1;
mtq_weight = 1;
thrust_weight = 1e-3;
rw_momentum_weight = 0.1;


%%

Yref_vec = [zeros(7,1);qRef(2:4);zeros(19,1)];
Yref_vec_tot = repmat(Yref_vec,prediction_horizon,1);

Yref_dotVec = [zeros(7,1);ones(3,1);zeros(19,1)];
Yref_dotVec_tot = repmat(Yref_dotVec,prediction_horizon,1);

u_mtq_dotVec = [ones(3,1);zeros(4,1)];
u_rw_dotVec = [zeros(3,1);ones(4,1)];
u_thrust_dotVec = ones(6,1);

u_mtq_dotVec_tot = repmat(u_mtq_dotVec,prediction_horizon,1);
u_rw_dotVec_tot = repmat(u_rw_dotVec,prediction_horizon,1);
u_thrust_dotVec_tot = repmat(u_thrust_dotVec,prediction_horizon,1);

U_lb_cont_tot = repmat(U_lb(1:7),prediction_horizon,1);
U_ub_cont_tot = repmat(U_ub(1:7),prediction_horizon,1);
U_lb_int_tot = repmat(U_lb(8:13),prediction_horizon,1);
U_ub_int_tot = repmat(U_ub(8:13),prediction_horizon,1);

Y_lb_tot = repmat(Y_lb,prediction_horizon,1);
Y_ub_tot = repmat(Y_ub,prediction_horizon,1);

Y_lb_dotVec_tot = repmat(Y_lb_dotVec,prediction_horizon,1);
Y_ub_dotVec_tot = repmat(Y_ub_dotVec,prediction_horizon,1);

Y_rw_momentum_dotVec = [zeros(13,1);ones(4,1);zeros(12,1)];
Y_rw_momentum_dotVec_tot = repmat(Y_rw_momentum_dotVec,prediction_horizon,1);

u_cont_integrated_mat = prediction_horizon.*eye(7);
u_int_integrated_mat = prediction_horizon.*eye(6);
for horizonIter = 1:(prediction_horizon-1)
    u_cont_integrated_mat = [ ...
        u_cont_integrated_mat, (prediction_horizon-horizonIter) .* eye(7)];
    u_int_integrated_mat = [ ...
        u_int_integrated_mat, (prediction_horizon-horizonIter) .* eye(6)];
end

%%

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

%%

prob = optimproblem('ObjectiveSense','min');

u_delta_int = optimvar('u_int',6*prediction_horizon,1,...
    'LowerBound',repmat(zeros(6,1),prediction_horizon,1), ...
    'UpperBound',U_ub_int_tot);
v_int = optimvar('v_int',6*prediction_horizon,1,'Type','integer','LowerBound',repmat(zeros(6,1),prediction_horizon,1),'UpperBound',repmat(ones(6,1),prediction_horizon,1));

u_delta_cont = optimvar('u_cont',7*prediction_horizon,1,...
    'LowerBound',U_lb_cont_tot,'UpperBound',U_ub_cont_tot);

Y_opt = optimvar('Y_opt',length(Yref_dotVec_tot),1);

%%


prob.Objective = 2*attitude_weight*Yref_dotVec_tot'* ...
        (temp_Y .* Y_opt) ...
        + ...
        2*mtq_weight* ones(1,7) * ...
        u_cont_integrated_mat * (temp_u_cont .* u_delta_cont .* u_mtq_dotVec_tot) ...
        + ...
        2*thrust_weight * ones(1,6) * ...
        u_int_integrated_mat * (temp_u_int .* u_delta_int) ...
        + ...
        2*rw_actuation_weight* ones(1,7)* ...
        u_cont_integrated_mat * (temp_u_cont .* u_delta_cont .* u_rw_dotVec_tot) ...
        + ...
        2*rw_momentum_weight*Y_rw_momentum_dotVec_tot'* ...
        (temp_Y .* Y_opt);
    
%         2*attitude_weight*Yref_dotVec_tot'* ...
%         (Y_opt - temp_Y); ...
%         + ...
%         2*mtq_weight* ones(1,7) * ...
%         u_cont_integrated_mat * (u_delta_cont .* u_mtq_dotVec_tot - ...
%         temp_u_cont .* u_mtq_dotVec_tot) ...
%         + ...
%         2*thrust_weight * ones(1,6) * ...
%         u_int_integrated_mat *  (u_delta_int - temp_u_int) ...
%         + ...
%         2*rw_actuation_weight* ones(1,7)* ...
%         u_cont_integrated_mat * (u_delta_cont .* u_rw_dotVec_tot - ...
%         temp_u_cont .* u_rw_dotVec_tot) ...
%         + ...
%         2*rw_momentum_weight*Y_rw_momentum_dotVec_tot'* ...
%         (Y_opt - temp_Y);

cons1 = Y_opt == chi_0 + B_chi_int * u_delta_int + B_chi_cont * u_delta_cont;
prob.Constraints.cons1 = cons1;

cons2 = Y_opt .* Y_lb_dotVec_tot >= Y_lb_tot;
prob.Constraints.cons2 = cons2;

cons3 = Y_opt .* Y_ub_dotVec_tot <= Y_ub_tot;
prob.Constraints.cons3 = cons3;

cons4 = u_delta_cont <= U_ub_cont_tot;
prob.Constraints.cons4 = cons4;

cons5 = u_delta_cont >= U_lb_cont_tot;
prob.Constraints.cons5 = cons5;

% cons6 = u_delta_cont <= temp_u_cont + delta_u_cont;
% prob.Constraints.cons6 = cons6;
% 
% cons7 = u_delta_cont >= temp_u_cont - delta_u_cont;
% prob.Constraints.cons7 = cons7;

prob.Constraints.u_int_maxconstr = u_delta_int <= U_ub_int_tot.*v_int;
prob.Constraints.u_int_minconstr = u_delta_int >= U_lb_int_tot.*v_int;

%show(prob);

options = optimoptions(@intlinprog, 'Display','off');
options = optimoptions(options,'LPOptimalityTolerance',1e-10,'RelativeGapTolerance',1e-8,...
                      'ConstraintTolerance',1e-9,'IntegerTolerance',1e-6); 
                  
                  
%%

prev_fval = 0;
[sol, fval]  = solve(prob,'options',options);

%thediff = 1e-9;
max_iterations = 100;
iter = 0; % iteration counter
                  
while iter < max_iterations
    % abs((fval - prev_fval)/prev_fval) > thediff

    iter = iter + 1;
    % fprintf(' Iteration: %i. \n', iter)

    prev_fval = fval;
    % Solve the problem with the new constraints
    [sol, fval] = solve(prob,'options',options);
%     sol.u_cont
%     sol.u_int
    
    gamma = 2/(iter + 2);
    temp_Y = temp_Y + ( sol.Y_opt - temp_Y )*gamma;
    temp_u_int = round( temp_u_int + ( sol.u_int - temp_u_int)*gamma ); 
    temp_u_cont = temp_u_cont + ( sol.u_cont - temp_u_cont )*gamma; 
    
    prob.Objective = 2*attitude_weight*Yref_dotVec_tot'* ...
        (temp_Y .* Y_opt) ...
        + ...
        2*mtq_weight* ones(1,7) * ...
        u_cont_integrated_mat * (temp_u_cont .* u_delta_cont .* u_mtq_dotVec_tot) ...
        + ...
        2*thrust_weight * ones(1,6) * ...
        u_int_integrated_mat * (temp_u_int .* u_delta_int) ...
        + ...
        2*rw_actuation_weight* ones(1,7)* ...
        u_cont_integrated_mat * (temp_u_cont .* u_delta_cont .* u_rw_dotVec_tot) ...
        + ...
        2*rw_momentum_weight*Y_rw_momentum_dotVec_tot'* ...
        (temp_Y .* Y_opt);
  
%     if iter>100
%         disp('here')
%     end
end

disp(['Number of iterations used: ',num2str(iter)])

cost_attitude = attitude_weight * dot( ...
        ( temp_Y .* Yref_dotVec_tot - Yref_vec_tot ).^2 ...
        , ...
        ones( size( Yref_vec_tot ) ) );
cost_mtq = mtq_weight * dot(u_cont_integrated_mat*temp_u_cont.^2, u_mtq_dotVec);
cost_rw_actuation = rw_actuation_weight * dot(u_cont_integrated_mat*temp_u_cont.^2, u_rw_dotVec);
cost_thrust = thrust_weight * ones(1,6) * u_int_integrated_mat*temp_u_int;
cost_actuation = cost_mtq + cost_rw_actuation + cost_thrust;
cost_rw_momentum = rw_momentum_weight * dot( (temp_Y .* Y_rw_momentum_dotVec_tot).^2, ...
    ones( size( Yref_vec_tot ) ) );
    
plotData.cost_attitude = [plotData.cost_attitude, cost_attitude];
plotData.cost_actuation = [plotData.cost_actuation, cost_actuation];
plotData.cost_rw_momentum = [plotData.cost_rw_momentum, cost_rw_momentum];

U = [temp_u_cont(1:7); temp_u_int(1:6)];


end

