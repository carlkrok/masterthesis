function [ eph ] = SimulateSatellite_RWPD( mjd0, satelliteFilename, ...
    timestep, duration, prediction_horizon, numControlVariables, ...
    numThrusters, numPropellant )


global MU_EARTH
global R_EARTH
global J2_EARTH

global simConfig
global rwData
global mtqData
global propulsionData
global satData
global missionData
global plotData
global satelliteConfiguration

run(satelliteFilename);

if simConfig.enableQuatRef
    K_p = 1e1;
    K_d = 1e1; % 10
else
    K_p = 0;
end

if simConfig.enableOmegaRef
    K_d = 1e1; % 1e1; % 
    K_dd = 1e-1; % 1e-1; % 
else
    K_dd = 0;
end


[r_PQW, v_PQW] = OrbitalElementsToPQW( MU_EARTH, sat.initCond.orb_h, ...
    sat.initCond.orb_e, sat.initCond.orb_T );
R_PQWToECI = RotMat_PQWToECI( sat.initCond.orb_i, sat.initCond.orb_O, ...
    sat.initCond.orb_w );
r_ECI = R_PQWToECI * r_PQW;
v_ECI = R_PQWToECI * v_PQW;

Y0 = [ r_ECI; v_ECI; sat.initCond.q_ECI; sat.initCond.w_body; sat.rw.w; ...
    sat.constr.Ixx_struct; sat.constr.Iyy_struct; sat.constr.Izz_struct; ...
    sat.constr.Ixy_struct; sat.constr.Ixz_struct; sat.constr.Iyz_struct]; 
for propIter = 1:numPropellant
    Y0 = [Y0; sat.propulsion.thrusters(propIter).rho]; 
end
Y0 = [Y0; sat.constr.body_init_com_struct];

    
U0 = zeros(numControlVariables,1);


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
    rwData.efficiency = sat.rw.efficiency;
    rwData.idlePower = sat.rw.idlePower;
else
    rwData.maxAcc = 0;
    rwData.maxVel = sat.rw.maxVel;
    rwData.maxVelocity = 0;
    rwData.efficiency = 1;
    rwData.idlePower = 0;
end

if simConfig.enableMTQ
    if not( sat.mtq.exists )
        error("Satellite not configured with Magnetorquer")
    end
    mtqData.maxDipoleMoment = sat.mtq.maxDipoleMoment;
    mtqData.powerFactor = sat.mtq.powerFactor;
else
    mtqData.maxDipoleMoment = 0;
    mtqData.powerFactor = Inf;
end

if simConfig.enablePropulsion
    if not( sat.propulsion.exists )
        error("Satellite not configured with Propulsion system")
    end
    propulsionData.maxThrust = sat.propulsion.maxThrust;
    propulsionData.power = sat.propulsion.power;
else
    propulsionData.minThrust = 0;
    propulsionData.maxThrust = 0;
    propulsionData.power = 0;
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

unitQuatTol = 1e-6;

eph = [0, Y0'];
Y = Y0';
U = U0;

omega_prev = Y0(11:13);

for currStep = 1:(round(duration/timestep))
    
    disp(['Iteration ',num2str(currStep),' of ',num2str(round(duration/timestep))]);

    thisT = currStep*timestep;
    mjd = missionData.mjd0 + thisT/86400;

    r_ECI = Y(1:3)';
    q_ECI = Y(7:10)';
    omega_body = Y(11:13)';
    rw_omega = Y(14:17)';
    I_mat_body = [Y(18), Y(21), Y(22); ...
                Y(21), Y(19), Y(23); ...
                Y(22), Y(23), Y(20)];

    I_mat_sys = I_mat_body + rwData.A_mat * rwData.I_mat * rwData.A_MPinv_mat;

    % Error variables describe rotation from current attitude to reference
    % attitude.
    
    if simConfig.enableQuatRef
        qError_ECI = QuaternionError(simConfig.referenceQuaternion, q_ECI);
    else
        qError_ECI = [1;0;0;0];
    end
    
    if simConfig.enableOmegaRef
        omegaError = omega_body - simConfig.referenceOmega;
        alphaError = omega_body - omega_prev;
        omega_prev = omega_body;
    else
        omegaError = omega_body;
        alphaError = zeros(3,1);
    end
    
    domega_rw = RW_NonlinPD( qError_ECI, omegaError, I_mat_sys, ...
    rw_omega, rwData.A_MPinv_mat, rwData.I_mat, rwData.maxAcc, ...
    rwData.maxVel, zeros(3,1), K_p, K_d, K_dd, alphaError );

 
    U = [zeros(3,1); domega_rw; zeros(numThrusters,1)];

    for uIter = 4:7
        if U(uIter) > rwData.maxAcc
            disp('Something is wrong! U lim exceeded')
            U(uIter) = rwData.maxAcc;
        elseif U(uIter) < -rwData.maxAcc
            disp('Something is wrong! U lim exceeded')
            U(uIter) = -rwData.maxAcc;
        end
    end
    
    plotData.mtq_m = [plotData.mtq_m; thisT, U(1:3)'];
    plotData.prop_f = [plotData.prop_f; thisT, U(8:8+numThrusters-1)'];
    
    % Implement first optimal control move
    opts = odeset('RelTol',1e-5,'AbsTol',1e-6);
    [t, thisEphY] = ode45( @(t, Y) SatelliteAcceleration(Y, U, ...
        timestep, mjd, numThrusters, numPropellant), [0; timestep], Y, opts );
    Y = thisEphY(length(thisEphY(:,1)),:);
    
    % Save plant states
    thisNorm = norm(Y(7:10));
    if abs(thisNorm - 1) > unitQuatTol
        %fprintf('Non-unit quaternion found in step %i of ephemeris, length %d - will be normalized.\n',stepIter,thisNorm )
        Y(7:10) = quatnormalize( Y(7:10) );
    end
    
    eph = [eph; thisT, Y];

end


end

