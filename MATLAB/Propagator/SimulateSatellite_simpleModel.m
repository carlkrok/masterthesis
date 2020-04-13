function [ eph ] = SimulateSatellite_simpleModel( satelliteFilename, mjd0, ...
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

if not( sat.rw.exists )
    error("Satellite not configured with Reaction Wheels")
end
rwData.A_mat = sat.rw.A_mat;
rwData.A_MPinv_mat = sat.rw.A_MPinv_mat;
rwData.I_mat = sat.rw.I_mat;
rwData.maxAcc = sat.rw.maxAcc;
rwData.maxVel = sat.rw.maxVel;
rwData.maxVelocity = sat.rw.maxVel;

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

nStates = 23;
nOutput = nStates;
nInputs = 13;
timeStep = 0.1;

SatNLMPC = nlmpc(nStates, nOutput, nInputs);
SatNLMPC.Ts = timeStep;
SatNLMPC.PredictionHorizon = 20;
SatNLMPC.ControlHorizon = 5;

SatNLMPC.Model.IsContinuousTime = true;
SatNLMPC.Model.NumberOfParameters = 4;
SatNLMPC.Model.StateFcn = @(Y, U, I_mat_body, com_struct, tot_mass, ...
    b_earth_eci) SatelliteAcceleration_simple(Y, U, ...
    I_mat_body, com_struct, tot_mass, b_earth_eci);

SatNLMPC.Optimization.CustomCostFcn = @(Y,U,e,data,I_mat_body, ...
    com_struct, tot_mass, b_earth_eci) SatCostFnc_simple(...
    Y,U,e,data,I_mat_body, com_struct, tot_mass, b_earth_eci);
SatNLMPC.Optimization.ReplaceStandardCost = true;
SatNLMPC.Optimization.UseSuboptimalSolution = true;

SatNLMPC.Optimization.SolverOptions.Algorithm = 'sqp';
SatNLMPC.Optimization.SolverOptions.SpecifyObjectiveGradient = true;
SatNLMPC.Optimization.SolverOptions.SpecifyConstraintGradient = true;
SatNLMPC.Optimization.SolverOptions.Display = 'final-detailed'; % 'iter-detailed'; %
SatNLMPC.Optimization.SolverOptions.MaxIter = 100;
SatNLMPC.Optimization.SolverOptions.OptimalityTolerance = 1.00e-06;
SatNLMPC.Optimization.SolverOptions.ConstraintTolerance = 1.00e-06;

for mtqVal = 1:3
    SatNLMPC.MV(mtqVal).Min = -satData.mtq.maxDipoleMoment;
    SatNLMPC.MV(mtqVal).Max = satData.mtq.maxDipoleMoment;
end

for rwVal = 4:7
    SatNLMPC.MV(rwVal).Min = -rwData.maxAcc;
    SatNLMPC.MV(rwVal).Max = rwData.maxAcc;
end

for rwVal = 14:17
    SatNLMPC.States(rwVal).Min = -rwData.maxVel;
    SatNLMPC.States(rwVal).Max = rwData.maxVel;
end

for propVal = 8:13
    SatNLMPC.MV(propVal).Min = 0; %satData.propulsion.minThrust;
    SatNLMPC.MV(propVal).Max = satData.propulsion.maxThrust;
end

for propVal = 18:23
    SatNLMPC.States(propVal).Min = 0; 
end

r0_ECI = Y0(1:3);
I_mat_body_0 = [Y0(18), Y0(21), Y0(22); ...
            Y0(21), Y0(19), Y0(23); ...
            Y0(22), Y0(23), Y0(20)];
com_struct_0 = Y0(30:32);
rho_thrusters_0 = Y0(24:29);

tot_mass_0 = satData.constr.body_mass;
for thrusterIter = 1:length(satData.propulsion.thrusters)
    [ thisCom, thisMass ] = propMassCoM( ...
        satData.propulsion.thrusters(thrusterIter).structureDim, ...
        rho_thrusters_0(thrusterIter) );
    tot_mass_0 = tot_mass_0 + thisMass;
end

r0_ECEF = posECIToECEF(0, missionData.mjd0, r0_ECI);
[ latitude_0, longitude_0, altitude_0 ] = ECEFToLLA( r0_ECEF );

NED_zAxis_0 = - r0_ECI./norm(r0_ECI);
NED_yAxis_0 = cross(NED_zAxis_0,[0;0;1]);
NED_xAxis_0 = cross(NED_yAxis_0, NED_zAxis_0);
rotMat_NEDToECI_0 = [ NED_xAxis_0, NED_yAxis_0, NED_zAxis_0 ];

b_earth_NED_0 = MagneticField( missionData.mjd0, latitude_0, longitude_0, altitude_0 );
b_earth_eci_0 = rotMat_NEDToECI_0*b_earth_NED_0;

Y0_simple = [Y0(1:17);Y0(24:29)];
validateFcns(SatNLMPC, Y0_simple, U0, [], { I_mat_body_0, com_struct_0,...
    tot_mass_0, b_earth_eci_0 });


SatNLMPC_options = nlmpcmoveopt;

unitQuatTol = 1e-6;

Duration = 3;
eph = [0, Y0'];
Y = Y0';
U = U0;
for currStep = 1:(Duration/timeStep)
    
    disp(['Iteration ',num2str(currStep),' of ',num2str(Duration/timeStep)]);
    
    thisT = currStep*timeStep;
    mjd = missionData.mjd0 + thisT/86400;
    
    r_ECI = Y(1:3)';
    I_mat_body = [Y(18), Y(21), Y(22); ...
                Y(21), Y(19), Y(23); ...
                Y(22), Y(23), Y(20)];
    com_struct = Y(30:32)';
    rho_thrusters = Y(24:29)';
    
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
    
    b_earth_NED = MagneticField( mjd, latitude, longitude, altitude );
    b_earth_eci = rotMat_NEDToECI*b_earth_NED;
    
    SatNLMPC_options.Parameters = { I_mat_body, com_struct, tot_mass, ...
        b_earth_eci };
    
    % Compute optimal control moves
    Y_simple = [Y(1:17),Y(24:29)];
    
    [U,SatNLMPC_options,info] = nlmpcmove(SatNLMPC,Y_simple,U,[],[],SatNLMPC_options)
    
    % Implement first optimal control move
    opts = odeset('RelTol',1e-5,'AbsTol',1e-6);
    [t, thisY] = ode45( @(t, Y) SatelliteAcceleration(Y, U, timeStep, mjd), [0; timeStep], Y, opts );
    Y = thisY(length(thisY(:,1)),:);
    
    % Save plant states
    thisNorm = norm(Y(7:10));
    if abs(thisNorm - 1) > unitQuatTol
        %fprintf('Non-unit quaternion found in step %i of ephemeris, length %d - will be normalized.\n',stepIter,thisNorm )
        Y(7:10) = quatnormalize( Y(7:10) );
    end
    
    eph = [eph; thisT, Y];
end



% prevEphTime = ephTimes(1);
% prevControllerTime = ephTimes(1);
% for ephTimeIter = 1:length(ephTimes)
%     nextEphTime = ephTimes(ephTimeIter);
%     ephStepLength = nextEphTime - prevEphTime;
%     
%     
% end
% eph = Ephemeris( Y0, stepTimes );





end

