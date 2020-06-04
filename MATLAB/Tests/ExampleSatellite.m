clear all

global satelliteConfiguration
global simConfig

% Choose satellite configuration
satelliteConfiguration = 2;

simConfig.enableJ2 = true;
simConfig.enableDrag = true;
simConfig.enableSRP = true;
simConfig.enableGravityGradient = true;

simConfig.enableRW = true;
%simConfig.enableRW = false;

simConfig.enableMTQ = true;
%simConfig.enableMTQ = false;

%simConfig.enablePropulsion = true;
simConfig.enablePropulsion = false;
%simConfig.enablePropulsionInertia = true;
simConfig.enablePropulsionInertia = false;

eulerFirst = pi/3;
eulerSecond = pi/8;
eulerThird = -pi/5;

simConfig.enablePointing = false;
%simConfig.pointingTarget_LLA = [63.4184922, 10.4005655, 0];
%simConfig.pointingTarget_ECEF = LLAToECEF( simConfig.pointingTarget_LLA );
simConfig.enableQuatRef = true;
simConfig.firstReferenceQuaternion = [1; 0; 0; 0];
simConfig.secondReferenceQuaternionTime = 30;
simConfig.secondReferenceQuaternion = EulerToQuaternion(eulerFirst,eulerSecond,eulerThird);
simConfig.thirdReferenceQuaternionTime = 90;
simConfig.thirdReferenceQuaternion = [1; 0; 0; 0];
simConfig.enableOmegaRef = false;
simConfig.firstReferenceOmega = [0; 0; 0];
simConfig.secondReferenceOmegaTime = 30;
simConfig.secondReferenceOmega = [0.0125; 0; 0];
simConfig.thirdReferenceOmegaTime = 90;
simConfig.thirdReferenceOmega = [0; 0; 0];


if satelliteConfiguration == 1
    satelliteFilename = "Satellite001";
    numControlVariables = 13;
    numThrusters = 6;
    numPropellant = 6;
elseif satelliteConfiguration == 2
    satelliteFilename = "Satellite002";
    numControlVariables = 19;
    numThrusters = 12;
    numPropellant = 1;
end

run earthParameters;
run(satelliteFilename);

t0_UTC = [2005 10 31 12 00 00.00];
t0_MJD = mjuliandate(t0_UTC);

global simConfig

global plotData
plotData.mtq_m = zeros(1,4);
plotData.prop_f = zeros(1,1+numThrusters);
plotData.cost_attitude = 0;
plotData.cost_actuation = 0;
plotData.cost_rw_momentum = 0;
plotData.disturbance_torques_norm = zeros(1,2);
plotData.srp_torque_norm = zeros(1,2);
plotData.drag_torque_norm = zeros(1,2);
plotData.grav_torque_norm = zeros(1,2);
plotData.wdot_sat_body = zeros(1,4);

sat_period = SatellitePeriod( MU_EARTH, sat.initCond.orb_a );
% stepLength = 1; % [s]
% numSteps = floor((0.05 * sat_period) / stepLength );
% stepTimes = zeros( 1, numSteps+1 );
% for stepIter = 0 : numSteps
%    stepTimes( stepIter+1 ) = (stepIter * stepLength);
% end

% MPC MTQ and nanofeep 5s * 20 = 100s pred
% MPC WaterJet 1s * 10 = 10s pred

timestep_pd = 0.1;
timestep_controller = 1; % 2; 
timestep_prediction = 1; % 2; 
prediction_horizon = 5;% 10; % 10 
duration = 150; %numSteps*stepLength;
eph = SimulateSatellite_integerMPC( t0_MJD, satelliteFilename, timestep_controller, ...
    timestep_prediction, duration, prediction_horizon, numControlVariables, ...
    numThrusters, numPropellant );
% eph = SimulateSatellite_RWPD( t0_MJD, satelliteFilename, ...
%     timestep_pd, duration, prediction_horizon, numControlVariables, ...
%     numThrusters, numPropellant );


timeVec = eph(:,1);

%%


% includeEarth = false;
% xyz = eph(:,2:4);
% PlotOrbit( xyz, includeEarth );
% 
% 
% %%
% 
% figure
% hold on
% grid on
% plot(1:length(omega_error_norm), omega_error_norm,'r')
% plot(1:length(omega_ref_norm), omega_ref_norm,'g')
% plot(1:length(omega_actual_norm), omega_actual_norm,'b')
% hold off
% 
% figure
% hold on
% grid on
% plot(1:length(omega_ref_norm), omega_ref(1,:),'r--')
% plot(1:length(omega_actual_norm), omega_actual(1,:),'r')
% plot(1:length(omega_ref_norm), omega_ref(2,:),'g--')
% plot(1:length(omega_actual_norm), omega_actual(2,:),'g')
% plot(1:length(omega_ref_norm), omega_ref(3,:),'b--')
% plot(1:length(omega_actual_norm), omega_actual(3,:),'b')
% hold off


% isCircular = true;
% xyz = eph(:,2:4);
% t = eph(:,1);
% PlotOrbitRadiusDeviation( t, xyz, isCircular );

   
%%

refTime = [0:simConfig.secondReferenceQuaternionTime, ...
    simConfig.secondReferenceQuaternionTime:simConfig.thirdReferenceQuaternionTime, ...
    simConfig.thirdReferenceQuaternionTime:duration];

xRef = [zeros(1,simConfig.secondReferenceQuaternionTime+1), ...
    eulerFirst .* ones(1,simConfig.thirdReferenceQuaternionTime-simConfig.secondReferenceQuaternionTime+1), ...
    zeros(1,duration-simConfig.thirdReferenceQuaternionTime+1)];

yRef = [zeros(1,simConfig.secondReferenceQuaternionTime+1), ...
    eulerSecond .* ones(1,simConfig.thirdReferenceQuaternionTime-simConfig.secondReferenceQuaternionTime+1), ...
    zeros(1,duration-simConfig.thirdReferenceQuaternionTime+1)];

zRef = [zeros(1,simConfig.secondReferenceQuaternionTime+1), ...
    eulerThird .* ones(1,simConfig.thirdReferenceQuaternionTime-simConfig.secondReferenceQuaternionTime+1), ...
    zeros(1,duration-simConfig.thirdReferenceQuaternionTime+1)];

quaternions = eph(:,8:11);
PlotEulerAngles( quaternions, timeVec, refTime, xRef, yRef, zRef );

quatRef = [simConfig.firstReferenceQuaternion * ones(1,simConfig.secondReferenceQuaternionTime+1), ...
    simConfig.secondReferenceQuaternion * ones(1,simConfig.thirdReferenceQuaternionTime-simConfig.secondReferenceQuaternionTime), ...
    simConfig.thirdReferenceQuaternion * ones(1,duration-simConfig.thirdReferenceQuaternionTime)];

PlotQuaternionError( quatRef, quaternions, timeVec)


%%

%PlotMultipleEulerAngles( eph3, eph2, eph )


%%

if simConfig.enableRW
rw_w = eph(:,15:18);
PlotRWMomentum( rw_w, timeVec )
end

%%


omega = eph(:,12:14);
PlotRotRate( omega, timeVec )


%%

% qArray = quaternion(eph(:,8:11));
% tp = theaterPlot;
% op = orientationPlotter(tp);
% axis([-1 1 -1 1 -1 1])
% view(3)
% for ephIter = 1:length(qArray)
%     %clf reset
%     plotOrientation(op,qArray(ephIter,:));
%     pause(1e-1)
% end


%%

if simConfig.enablePropulsion
prop_rho = eph(:,25:25+numPropellant-1);
PlotPropRho( prop_rho, timeVec )
PlotPropMass( prop_rho, timeVec )
PlotTotPropMass( prop_rho, timeVec )
com = eph(:,25+numPropellant:25+numPropellant+2);
PlotCoM( com, timeVec )
end

%%


PlotPropForce( plotData.prop_f );

PlotMTQDipole( plotData.mtq_m )

% PlotCosts(plotData.cost_attitude,plotData.cost_actuation,plotData.cost_rw_momentum);

% PlotDisturbanceTorqueNorm ( plotData.disturbance_torques_norm, ...
%     plotData.srp_torque_norm, plotData.drag_torque_norm, ...
%     plotData.grav_torque_norm);

PlotOmegaDot( plotData.wdot_sat_body );

%%

%save('e1')
