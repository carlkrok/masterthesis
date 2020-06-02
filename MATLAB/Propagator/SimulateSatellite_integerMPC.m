function [ eph ] = SimulateSatellite_integerMPC( mjd0, satelliteFilename, timestep_controller, ...
    timestep_prediction, duration, prediction_horizon, numControlVariables, ...
    numThrusters, numPropellant)


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
rw_vel_ref = 1000*2*pi/60;

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
    propulsionData.thrustRange = sat.propulsion.thrustRange;
    propulsionData.power = sat.propulsion.power;
else
    propulsionData.minThrust = 0;
    propulsionData.thrustRange = 0;
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

U_lb = [ -mtqData.maxDipoleMoment .* ones(3,1);...
        -rwData.maxAcc .* ones(4,1); ...
        ones(numThrusters,1)];
U_ub = [ mtqData.maxDipoleMoment .* ones(3,1); ...
    rwData.maxAcc .* ones(4,1); ...
    ones(numThrusters,1).*length(propulsionData.thrustRange)]; % propulsionData.maxThrust .* 

Y_lb = [zeros(7,1); -rwData.maxVel .*ones(4,1); zeros(numPropellant,1)];
Y_lb_dotVec = [zeros(7,1); ones(4,1); zeros(numPropellant,1)];
    
Y_ub = [zeros(7,1); rwData.maxVel .*ones(4,1); zeros(numPropellant,1)];
Y_ub_dotVec = [zeros(7,1); ones(4,1); zeros(numPropellant,1)];

for currStep = 1:(round(duration/timestep_controller)) % -1) % DEBUG
    
    disp(['Iteration ',num2str(currStep),' of ',num2str(round(duration/timestep_controller))]);

    thisT = currStep*timestep_controller;
    mjd = missionData.mjd0 + thisT/86400;
        
%     if currStep == 1 % DEBUG


        r_ECI = Y(1:3)';
        q_ECI = Y(7:10)';
        omega_body = Y(11:13)';
        I_mat_body = [Y(18), Y(21), Y(22); ...
                    Y(21), Y(19), Y(23); ...
                    Y(22), Y(23), Y(20)];
        com_struct = Y(24+numPropellant:24+numPropellant+2)';
        rho_thrusters = Y(24:24+numPropellant-1)';

        I_mat_sys = I_mat_body + rwData.A_mat * rwData.I_mat * rwData.A_MPinv_mat;

        tot_mass = satData.constr.body_mass;
        for thrusterIter = 1:numPropellant
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
        
%         global B_EARTH_ECI_DEBUG
%         B_EARTH_ECI_DEBUG = b_earth_eci;
    
        if simConfig.enableRW
            B_RW = LinearizedRW( I_mat_sys, rwData.I_mat, rwData.A_mat );
        else
            B_RW = zeros(3,4);
        end

        if simConfig.enableMTQ
            B_MTQ = LinearizedMTQ( I_mat_sys, b_earth_body );
        else
            B_MTQ = zeros(3,3);
        end

        if simConfig.enablePropulsion
            B_PROP = LinearizedPROP( I_mat_sys, com_struct, ...
                satData.propulsion.thrusters, numThrusters );
        else
            B_PROP = zeros(3,numThrusters);
        end

        
        A_cont = LinearizedStateMat( q_ECI, numPropellant );
        B_cont = [zeros(4,8+numThrusters-1); ...
        B_MTQ, B_RW, B_PROP; ...
        zeros(4,3), eye(4,4), zeros(4,numThrusters); ...
        zeros(numPropellant,7), LinearizedRHO( ...
            satData.propulsion.thrusters, numThrusters)];

        A_disc = eye(size(A_cont)) + timestep_prediction .* A_cont;
        B_disc = timestep_prediction .* B_cont;



        Y_simple = [Y(7:17),Y(24:24+numPropellant-1)]';
        
        
        
        
        
        predT = thisT + timestep_prediction;
        
        qRef = simConfig.firstReferenceQuaternion;
        omegaRef = simConfig.firstReferenceOmega;

        if predT > simConfig.secondReferenceQuaternionTime
            qRef = simConfig.secondReferenceQuaternion;
        end
        if predT > simConfig.thirdReferenceQuaternionTime
            qRef = simConfig.thirdReferenceQuaternion;
        end

        if predT > simConfig.secondReferenceOmegaTime
            omegaRef = simConfig.secondReferenceOmega;
        end
        if predT > simConfig.thirdReferenceOmegaTime
            omegaRef = simConfig.thirdReferenceOmega;
        end

        chi_ref = [0;qRef(2:4);omegaRef;rw_vel_ref.*ones(4,1); ...
            zeros(numPropellant,1)];
            
        for chiIter = 2:prediction_horizon
            
            predT = thisT + timestep_prediction*chiIter;
            
            qRef = simConfig.firstReferenceQuaternion;
            omegaRef = simConfig.firstReferenceOmega;
            
            if predT > simConfig.secondReferenceQuaternionTime
                qRef = simConfig.secondReferenceQuaternion;
            end
            if predT > simConfig.thirdReferenceQuaternionTime
                qRef = simConfig.thirdReferenceQuaternion;
            end
            
            if predT > simConfig.secondReferenceOmegaTime
                omegaRef = simConfig.secondReferenceOmega;
            end
            if predT > simConfig.thirdReferenceOmegaTime
                omegaRef = simConfig.thirdReferenceOmega;
            end
        
            chi_ref = [chi_ref;[0;qRef(2:4);omegaRef;rw_vel_ref.*ones(4,1); ...
                zeros(numPropellant,1)]];
            
        end
        
        
    
%         Uall = ga_MPC(Y_simple, A_disc,B_disc,prediction_horizon, ...
%             numControlVariables, numThrusters, U_lb, U_ub, ...
%             Y_lb, Y_ub, Y_lb_dotVec, Y_ub_dotVec, chi_ref, rw_vel_ref);
%         
%     end
%     
%     U = Uall(1+currStep*13:13+currStep*13);
    
    

    U = ga_MPC(Y_simple, A_disc,B_disc,prediction_horizon, ...
        numControlVariables, numThrusters, numPropellant, U_lb, U_ub, ...
        Y_lb, Y_ub, Y_lb_dotVec, Y_ub_dotVec, chi_ref, rw_vel_ref);


    U(8:8+numThrusters-1) = propulsionData.thrustRange(U(8:8+numThrusters-1));
    
    
    for uIter = 1:3
        if U(uIter) > mtqData.maxDipoleMoment
            disp('Something is wrong! U lim exceeded')
            U(uIter) = mtqData.maxDipoleMoment;
        elseif U(uIter) < -mtqData.maxDipoleMoment
            disp('Something is wrong! U lim exceeded')
            U(uIter) = -mtqData.maxDipoleMoment;
        end
    end

    for uIter = 4:7
        if U(uIter) > rwData.maxAcc
            disp('Something is wrong! U lim exceeded')
            U(uIter) = rwData.maxAcc;
        elseif U(uIter) < -rwData.maxAcc
            disp('Something is wrong! U lim exceeded')
            U(uIter) = -rwData.maxAcc;
        end
    end

    for uIter = 8:8+numThrusters-1
        if U(uIter) > propulsionData.maxThrust
            disp('Something is wrong! U lim exceeded')
            U(uIter) = propulsionData.maxThrust;
        elseif U(uIter) < 0
            disp('Something is wrong! U lim exceeded')
            U(uIter) = 0;
        end
    end
    
    plotData.mtq_m = [plotData.mtq_m; thisT, U(1:3)'];
    plotData.prop_f = [plotData.prop_f; thisT, U(8:8+numThrusters-1)'];
    
    % Implement first optimal control move
    opts = odeset('RelTol',1e-5,'AbsTol',1e-6);
    [t, thisEphY] = ode45( @(t, Y) SatelliteAcceleration(Y, U, ...
        timestep_controller, mjd, numThrusters, numPropellant), [0; timestep_controller], Y, opts );
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

