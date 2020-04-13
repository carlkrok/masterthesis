function [ eph ] = SimulateSatellite_fullModel( satelliteFilename, mjd0, ...
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

nStates = 32;
nOutput = nStates;
nInputs = 13;
timeStep = 1;

SatNLMPC = nlmpc(nStates, nOutput, nInputs);
SatNLMPC.Ts = timeStep;
SatNLMPC.PredictionHorizon = 8;
SatNLMPC.ControlHorizon = 2;

SatNLMPC.Model.IsContinuousTime = true;
SatNLMPC.Model.NumberOfParameters = 2;
SatNLMPC.Model.StateFcn = @(Y, U, Ts, mjd) SatelliteAcceleration(Y, U, Ts, mjd);

SatNLMPC.Optimization.CustomCostFcn = @(Y,U,e,data,Ts, mjd) SatCostFnc(Y,U,e,data,Ts, mjd);
SatNLMPC.Optimization.ReplaceStandardCost = true;
SatNLMPC.Optimization.UseSuboptimalSolution = true;

SatNLMPC.Optimization.SolverOptions.Algorithm = 'sqp';
%SatNLMPC.Optimization.SolverOptions.SpecifyObjectiveGradient = true;
%SatNLMPC.Optimization.SolverOptions.SpecifyConstraintGradient = true;
SatNLMPC.Optimization.SolverOptions.Display = 'final-detailed'; % 'iter-detailed'; %
SatNLMPC.Optimization.SolverOptions.MaxIter = 5;
SatNLMPC.Optimization.SolverOptions.OptimalityTolerance = 1.00e-05;
SatNLMPC.Optimization.SolverOptions.ConstraintTolerance = 1.00e-05;

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

for propVal = 24:29
    SatNLMPC.States(propVal).Min = 0; 
end



validateFcns(SatNLMPC, Y0, U0, [], {timeStep, missionData.mjd0});


SatNLMPC_options = nlmpcmoveopt;

unitQuatTol = 1e-6;

Duration = 200;
eph = [0, Y0'];
Y = Y0;
U = U0;
for currStep = 1:(Duration/timeStep)
    thisT = currStep*timeStep;
    mjd = missionData.mjd0 + thisT/86400;
    SatNLMPC_options.Parameters = {thisT, mjd};
    % Compute optimal control moves
    [U,SatNLMPC_options, info] = nlmpcmove(SatNLMPC,Y,U,[],[],SatNLMPC_options)
    % Implement first optimal control move
    opts = odeset('RelTol',1e-5,'AbsTol',1e-6,'OutputFcn',@odeprog,'Events',@odeabort);
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

