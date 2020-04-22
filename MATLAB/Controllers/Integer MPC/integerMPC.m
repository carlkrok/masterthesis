function U = integerMPC(Y0, A_disc,B_disc,prediction_horizon, U_lb, U_ub, ...
    Y_lb, Y_ub, Y_lb_dotVec, Y_ub_dotVec, qRef)

global plotData

temp_Y = repmat(Y0, prediction_horizon, 1);
temp_u_int = repmat(zeros(6,1), prediction_horizon, 1);
temp_u_cont = repmat(zeros(7,1), prediction_horizon, 1);
delta_u_cont = repmat([0.1;0.1;0.1;1;1;1;1],prediction_horizon,1);

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

u_delta_int = optimvar('u_int',6*prediction_horizon,1,...
    'LowerBound',repmat(zeros(6,1),prediction_horizon,1), ...
    'UpperBound',U_ub_int_tot);
v_int = optimvar('v_int',6*prediction_horizon,1,'Type','integer','LowerBound',0,'UpperBound',1);

u_delta_cont = optimvar('u_cont',7*prediction_horizon,1,...
    'LowerBound',U_lb_cont_tot,'UpperBound',U_ub_cont_tot);

Y_opt = optimvar('Y_opt',length(Yref_dotVec_tot),1);

z_total = optimvar('z_total',1,'LowerBound',0);
% z_attitude = optimvar('z_attitude',1,'LowerBound',0);
% z_mtq = optimvar('z_mtq',1,'LowerBound',0);
% z_rw_actuation = optimvar('z_rw_actuation',1,'LowerBound',0);
% z_rw_momentum = optimvar('z_rw_momentum',1,'LowerBound',0);
% z_thrust = optimvar('z_thrust',1,'LowerBound',0);

u_cont_integrated_mat = prediction_horizon.*eye(7);
u_int_integrated_mat = prediction_horizon.*eye(6);
for horizonIter = 1:(prediction_horizon-1)
    u_cont_integrated_mat = [ ...
        u_cont_integrated_mat, (prediction_horizon-horizonIter) .* eye(7)];
    u_int_integrated_mat = [ ...
        u_int_integrated_mat, (prediction_horizon-horizonIter) .* eye(6)];
end

attitude_weight = 1e6;
rw_actuation_weight = 1;
mtq_weight = 1;
thrust_weight = 1e-3;
rw_momentum_weight = 0.1;

prob.Objective = z_total;
%         attitude_weight * z_attitude; ...
%     + mtq_weight * z_mtq ...
%     + rw_actuation_weight * z_rw_actuation ...
%     + thrust_weight * ones(1,6) * u_int_integrated_mat * u_delta_int ... * z_thrust ...
%     + rw_momentum_weight * z_rw_momentum;

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

cons6 = u_delta_cont <= temp_u_cont + delta_u_cont;
prob.Constraints.cons6 = cons6;

cons7 = u_delta_cont >= temp_u_cont - delta_u_cont;
prob.Constraints.cons7 = cons7;

prob.Constraints.u_int_maxconstr = u_delta_int <= U_ub_int_tot.*v_int;
prob.Constraints.u_int_minconstr = u_delta_int >= U_ub_int_tot.*v_int;

%show(prob);

options = optimoptions(@intlinprog, 'Display','off');
options = optimoptions(options,'LPOptimalityTolerance',1e-10,'RelativeGapTolerance',1e-8,...
                      'ConstraintTolerance',1e-9,'IntegerTolerance',1e-6); 
                  
sol  = solve(prob,'options',options);

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

cost_total = cost_attitude + cost_mtq + cost_rw_actuation + cost_thrust + ...
    cost_rw_momentum;

thediff = 1e-4;
iter = 1; % iteration counter
temp_Y = sol.Y_opt;
temp_u_int = sol.u_int;
temp_u_cont = sol.u_cont;

zslack_total = sol.z_total;
%zslack_attitude = sol.z_attitude;
%zslack_mtq = sol.z_mtq;
%zslack_rw_actuation = sol.z_rw_actuation;
%zslack_rw_momentum = sol.z_rw_momentum;
%zslack_thrust = sol.z_thrust;
                  
while abs((zslack_total - cost_attitude)/cost_attitude) > thediff
        %abs((zslack_attitude - cost_attitude)/cost_attitude) > thediff %|| ...
        %abs((zslack_mtq - cost_mtq)/cost_mtq) > thediff || ...
        %abs((zslack_rw_actuation - cost_rw_actuation)/cost_rw_actuation) > thediff || ...
        %abs((zslack_rw_momentum - cost_rw_momentum)/cost_rw_momentum) > thediff %|| ...
        %abs((zslack_thrust - cost_thrust)/cost_thrust) > thediff
    % relative error
    
    fprintf(' Iteration: %i. \n', iter)
%     fprintf(' Iteration: %i. \n Relative attitude error: %d. \n Relative mtq error: %d. \n Relative rw_actuation error: %d. \n Relative rw_moment error: %d. \n', ...
%         iter, abs((zslack_attitude - cost_attitude)/cost_attitude), ...
%         abs((zslack_mtq - cost_mtq)/cost_mtq), ...
%         abs((zslack_rw_actuation - cost_rw_actuation)/cost_rw_actuation), ...
%         abs((zslack_rw_momentum - cost_rw_momentum)/cost_rw_momentum))
%         %abs((zslack_thrust - cost_thrust)/cost_thrust) ... \n Relative thrust error: %d 
    
    constr_total = 2*attitude_weight*(temp_Y .* Yref_dotVec_tot)'* ...
        (Y_opt .* Yref_dotVec_tot) + ...
        - z_total <= cost_total;
%         2*mtq_weight* ones(1,7) * ...
%         u_cont_integrated_mat * (u_delta_cont .* temp_u_cont .* ...
%         u_mtq_dotVec_tot) + ...
%         2*thrust_weight * ones(1,6) * ...
%         u_int_integrated_mat * temp_u_int .* u_delta_int + ...
%         2*rw_actuation_weight* ones(1,7)* ...
%         u_cont_integrated_mat * (u_delta_cont .* temp_u_cont .* ...
%         u_rw_dotVec_tot) + ...
%         2*rw_momentum_weight*(temp_Y .* ...
%         Y_rw_momentum_dotVec_tot)'* ...
%         (Y_opt .* Y_rw_momentum_dotVec_tot) ...
            
    newname_constr_total = ['total_iteration',num2str(iter)];
    prob.Constraints.(newname_constr_total) = constr_total;
    
%     constr_attitude = 2*attitude_weight*(temp_Y .* Yref_dotVec_tot)'* ...
%         (Y_opt .* Yref_dotVec_tot) - z_attitude <= cost_attitude;
%     newname_constr_attitude = ['attitude_iteration',num2str(iter)];
%     prob.Constraints.(newname_constr_attitude) = constr_attitude;
    
%     constr_mtq = 2*mtq_weight* ones(1,7) * ...
%         u_cont_integrated_mat * (u_delta_cont .* temp_u_cont .* ...
%         u_mtq_dotVec_tot) - z_mtq <= cost_mtq;
%     newname_constr_mtq = ['mtq_iteration',num2str(iter)];
%     prob.Constraints.(newname_constr_mtq) = constr_mtq;
    
%     constr_thrust = 2*thrust_weight * ones(1,6) * ...
%         u_int_integrated_mat * temp_u_int .* u_delta_int - z_thrust ...
%         <= cost_thrust;
%     newname_constr_thrust = ['thrust_iteration',num2str(iter)];
%     prob.Constraints.(newname_constr_thrust) = constr_thrust;
    
%     constr_rw_actuation = 2*rw_actuation_weight* ones(1,7)* ...
%         u_cont_integrated_mat * (u_delta_cont .* temp_u_cont .* ...
%         u_rw_dotVec_tot) - z_rw_actuation <= cost_rw_actuation;
%     newname_constr_rw_actuation = ['rw_actuation_iteration',num2str(iter)];
%     prob.Constraints.(newname_constr_rw_actuation) = constr_rw_actuation;

%     constr_rw_momentum = 2*rw_momentum_weight*(temp_Y .* ...
%         Y_rw_momentum_dotVec_tot)'* ...
%         (Y_opt .* Y_rw_momentum_dotVec_tot) - z_rw_momentum ...
%         <= cost_rw_momentum;
%     newname_constr_rw_momentum = ['rw_momentum_iteration',num2str(iter)];
%     prob.Constraints.(newname_constr_rw_momentum) = constr_rw_momentum;

    % Solve the problem with the new constraints
    sol = solve(prob,'options',options);
    
    if any(isempty(sol))
        sol.Y_opt = tempY;
        sol.u_int = temp_u_int;
        sol.u_cont = temp_u_cont;
    end
    
%    assets = (assets+xLinInt.xvars)/2; % Midway from the previous to the current
%    assets = xLinInt(xvars); % Use the previous line or this one
    
    temp_Y =  sol.Y_opt; % (temp_Y + sol.Y_opt)/2; % 
    temp_u_int = sol.u_int; % ( temp_u_int + sol.u_int)/2; % 
    temp_u_cont = sol.u_cont; % ( temp_u_cont + sol.u_cont)/2; %
    
%     cons6 = u_delta_cont <= temp_u_cont + delta_u_cont;
%     prob.Constraints.cons6 = cons6;
% 
%     cons7 = u_delta_cont >= temp_u_cont - delta_u_cont;
%     prob.Constraints.cons7 = cons7;

    zslack_total = sol.z_total; % ( zslack_total + )/2
%    zslack_attitude = sol.z_attitude; %( zslack_attitude + sol.z_attitude)/2;
%    zslack_mtq = sol.z_mtq; %( zslack_mtq + sol.z_mtq)/2;
%    zslack_rw_actuation = sol.z_rw_actuation; %( zslack_rw_actuation + sol.z_rw_actuation)/2;
%    zslack_rw_momentum = sol.z_rw_momentum; %( zslack_rw_momentum + sol.z_rw_momentum)/2;
%    zslack_thrust = sol.z_thrust; %( zslack_thrust + sol.z_thrust)/2;
    
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
    
    cost_total = cost_attitude + cost_mtq + cost_rw_actuation + cost_thrust + ...
    cost_rw_momentum

    cons6 = u_delta_cont <= temp_u_cont + delta_u_cont;
    prob.Constraints.cons6 = cons6;

    cons7 = u_delta_cont >= temp_u_cont - delta_u_cont;
    prob.Constraints.cons7 = cons7;

    iter = iter + 1;
end

disp(['Number of iterations used: ',num2str(iter)])

plotData.cost_attitude = [plotData.cost_attitude, cost_attitude];
plotData.cost_actuation = [plotData.cost_actuation, cost_actuation];
plotData.cost_rw_momentum = [plotData.cost_rw_momentum, cost_rw_momentum];

U = [temp_u_cont(1:7); temp_u_int(1:6)];


end

