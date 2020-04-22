function [ eph ] = SimulateSatellite_linearizedCustom( satelliteFilename, mjd0, ...
    ephTimes, controllerFreq )

global MU_EARTH
global R_EARTH
global J2_EARTH

global simConfig
global rwData
global mtqData
global propulsionData
global satData
global missionData

global omega_error_norm
global omega_ref_norm
global omega_actual_norm
global omega_error
global omega_ref
global omega_actual

K_p = 0.1;
K_d = 0.3;

run(satelliteFilename);

[r_PQW, v_PQW] = OrbitalElementsToPQW( MU_EARTH, sat.initCond.orb_h, ...
    sat.initCond.orb_e, sat.initCond.orb_T );
R_PQWToECI = RotMat_PQWToECI( sat.initCond.orb_i, sat.initCond.orb_O, ...
    sat.initCond.orb_w );
r_ECI = R_PQWToECI * r_PQW;
v_ECI = R_PQWToECI * v_PQW;

Y0 = [ r_ECI; v_ECI; sat.initCond.q_ECI; sat.initCond.w_body; sat.rw.w; ...
    sat.constr.Ixx_struct; sat.constr.Iyy_struct; sat.constr.Izz_struct; ...
    sat.constr.Ixy_struct; sat.constr.Ixz_struct; sat.constr.Iyz_struct; ...
    sat.propulsion.thrusters(1).rho; ...
    sat.propulsion.thrusters(2).rho; ...
    sat.propulsion.thrusters(3).rho; ...
    sat.propulsion.thrusters(4).rho; ...
    sat.propulsion.thrusters(5).rho; ...
    sat.propulsion.thrusters(6).rho; ...
    sat.constr.body_init_com_struct];

U0 = zeros(13,1);


rwData.A_mat = sat.rw.A_mat;
rwData.A_MPinv_mat = sat.rw.A_MPinv_mat;
rwData.I_mat = sat.rw.I_mat;
if simConfig.enableRW
    if not( sat.rw.exists )
        error("Satellite not configured with Reaction Wheels")
    end
    rwData.maxAcc = sat.rw.maxAcc;
    rwData.maxVel = sat.rw.maxVel;
    rwData.maxVelocity = sat.rw.maxVel;
else
    rwData.maxAcc = 0;
    rwData.maxVel = 0;
    rwData.maxVelocity = 0;
end

if simConfig.enableMTQ
    if not( sat.mtq.exists )
        error("Satellite not configured with Magnetorquer")
    end
    mtqData.maxDipoleMoment = sat.mtq.maxDipoleMoment;
else
    mtqData.maxDipoleMoment = 0;
end

if simConfig.enablePropulsion
    if not( sat.propulsion.exists )
        error("Satellite not configured with Propulsion system")
    end
    propulsionData.minThrust = sat.propulsion.minThrust;
    propulsionData.maxThrust = sat.propulsion.maxThrust;
else
    propulsionData.minThrust = 0;
    propulsionData.maxThrust = 0;
end

satData.constr = sat.constr;
satData.propulsion = sat.propulsion;
satData.mtq = sat.mtq;
satData.surfaceCenterVectorsNormalVectorsAreas = sat.constr.surfaceCenterVectorsNormalVectorsAreas;
satData.coeffR = sat.constr.coeffR;
satData.coeffDrag = sat.constr.coeffDrag;

missionData.mjd0 = mjd0;
missionData.mu = MU_EARTH;
missionData.Re = R_EARTH;
missionData.J2 = J2_EARTH;


%%

timeStep_discretization = 0.1;
timeStep_sim = 0.1;

C = eye(17,17);
D = zeros(17,13);

%%

unitQuatTol = 1e-6;

Duration = 10;
eph = [0, Y0'];
Y = Y0';
U = U0;

for currStep = 1:(Duration/timeStep_sim)
    
    disp(['Iteration ',num2str(currStep),' of ',num2str(Duration/timeStep_sim)]);

    thisT = currStep*timeStep_sim;
    mjd = missionData.mjd0 + thisT/86400;
    
    r_ECI = Y(1:3)';
    q_ECI = Y(7:10)';
    omega_body = Y(11:13)';
    I_mat_body = [Y(18), Y(21), Y(22); ...
                Y(21), Y(19), Y(23); ...
                Y(22), Y(23), Y(20)];
    com_struct = Y(30:32)';
    rho_thrusters = Y(24:29)';
    
    I_mat_sys = I_mat_body + rwData.A_mat * rwData.I_mat * rwData.A_MPinv_mat;
    
    tot_mass = satData.constr.body_mass;
    for thrusterIter = 1:length(satData.propulsion.thrusters)
        [ thisCom, thisMass ] = propMassCoM( ...
            satData.propulsion.thrusters(thrusterIter).structureDim, ...
            rho_thrusters(thrusterIter) );
        tot_mass = tot_mass + thisMass;
    end
    
    r_ECEF = posECIToECEF(thisT, missionData.mjd0, r_ECI);
    [ latitude, longitude, altitude ] = ECEFToLLA( r_ECEF );

    NED_zAxis = - r_ECI./norm(r_ECI);
    NED_yAxis = cross(NED_zAxis,[0;0;1]);
    NED_xAxis = cross(NED_yAxis, NED_zAxis);
    rotMat_NEDToECI = [ NED_xAxis, NED_yAxis, NED_zAxis ];
    
    rotMat_ECIToBody = QuaternionToRotMat(q_ECI);
    
    b_earth_NED = MagneticField( mjd, latitude, longitude, altitude );
    b_earth_eci = rotMat_NEDToECI*b_earth_NED;
    b_earth_body = rotMat_ECIToBody * b_earth_eci;
    
    A = LinearizedStateMat( r_ECI, MU_EARTH, q_ECI );
    B = [zeros(10,13); ...
    LinearizedMTQ( I_mat_sys, b_earth_body ), ...
    LinearizedRW( I_mat_sys, rwData.I_mat, rwData.A_mat ), ...
    LinearizedPROP( I_mat_sys, com_struct, satData.propulsion.thrusters ); ...
    zeros(4,3), eye(4,4), zeros(4,6) ];

    plant_cont = ss(A,B,C,D);
    plant_disc = c2d(plant_cont,timeStep_discretization);

    qError_ECI = QuaternionError(simConfig.referenceQuaternion', q_ECI);

    omega_ref_sat = - K_p * qError_ECI(2:4) - K_d * omega_body;
    
    omega_ref_controller = omega_ref_sat - omega_body;
    
    lb = [ -satData.mtq.maxDipoleMoment .* ones(3,1);...
        -rwData.maxAcc .* ones(4,1); zeros(6,1); 0];
    ub = [ satData.mtq.maxDipoleMoment .* ones(3,1); ...
        rwData.maxAcc .* ones(4,1); ...
        satData.propulsion.maxThrust .* ones(6,1); 1];
    
%     U = lsqlin(plant_disc.B(11:13,:),omega_ref_controller,[],[],[],[],lb,ub);

%     ControlDistMat = plant_disc.B(11:13,:)'/(plant_disc.B(11:13,:)*plant_disc.B(11:13,:)');
%     U = ControlDistMat * omega_ref_controller;

    prob = optimproblem('ObjectiveSense','min');
    
    u = optimvar('u',14,1,'LowerBound',lb,'UpperBound',ub);
    
    prob.Objective = dot((plant_disc.B(11:13,:)*u(1:13) - omega_ref_controller).^2,[1,1,1]);
    
    cons1 = plant_disc.B(11:13,:)*u(1:13) == u(14)*omega_ref_controller;
    prob.Constraints.cons1 = cons1;
    
    %show(prob)
    
    sol = solve(prob);
    
    U = sol.u(1:13);
    
    resulting_omega = plant_disc.B(11:13,:)*U;
    
    omega_error_norm = [omega_error_norm, norm(omega_ref_controller-resulting_omega)];
    omega_ref_norm = [omega_ref_norm, norm(omega_ref_controller)];
    omega_actual_norm = [omega_actual_norm, norm(resulting_omega)];
    
    omega_error = [omega_error, norm(omega_ref_controller-resulting_omega)];
    omega_ref = [omega_ref, omega_ref_controller];
    omega_actual = [omega_actual, resulting_omega];
    
    
    % Implement first optimal control move
    opts = odeset('RelTol',1e-5,'AbsTol',1e-6);
    [t, thisY] = ode45( @(t, Y) SatelliteAcceleration(Y, U, timeStep_sim, mjd), [0; timeStep_sim], Y, opts );
    Y = thisY(length(thisY(:,1)),:);
    
    % Save plant states
    thisNorm = norm(Y(7:10));
    if abs(thisNorm - 1) > unitQuatTol
        %fprintf('Non-unit quaternion found in step %i of ephemeris, length %d - will be normalized.\n',stepIter,thisNorm )
        Y(7:10) = quatnormalize( Y(7:10) );
    end
    
    eph = [eph; thisT, Y];
end


end

