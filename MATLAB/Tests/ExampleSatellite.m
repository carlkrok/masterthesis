clear all

satelliteFilename = "Satellite001";

run earthParameters;
run(satelliteFilename);

global simConfig

t0_UTC = [2005 10 31 12 56 37.67];
t0_MJD = mjuliandate(t0_UTC);

sat_period = SatellitePeriod( MU_EARTH, sat.initCond.orb_a );
stepLength = 1; % [s]
numSteps = floor((0.05 * sat_period) / stepLength );
stepTimes = zeros( 1, numSteps+1 );
for stepIter = 0 : numSteps
   stepTimes( stepIter+1 ) = (stepIter * stepLength);
end

simConfig.enableJ2 = true;
simConfig.enableDrag = true;
simConfig.enableSRP = true;
simConfig.enableGravityGradient = true;

simConfig.enableRW = true;
%simConfig.enableRW = false;

simConfig.enableMTQ = true;
%simConfig.enableMTQ = false;

simConfig.enablePropulsion = true;
%simConfig.enablePropulsion = false;

simConfig.enablePointing = false;
%simConfig.pointingTarget_LLA = [63.4184922, 10.4005655, 0];
%simConfig.pointingTarget_ECEF = LLAToECEF( simConfig.pointingTarget_LLA );
simConfig.referenceQuaternion = [1; 0; 0; 0];

%eph = SimulateSatellite_fullModel( satelliteFilename, t0_MJD, stepTimes );
%eph = SimulateSatellite_simpleModel( satelliteFilename, t0_MJD, stepTimes );
%eph = SimulateSatellite_linearizedMatlab( satelliteFilename, t0_MJD, stepTimes );
%eph = SimulateSatellite_linearizedCustom( satelliteFilename, t0_MJD, stepTimes );

timestep = 1;
prediction_horizon = 5;
duration = 50; %numSteps*stepLength;
eph = SimulateSatellite_customMPC( satelliteFilename, t0_MJD, ...
    timestep, duration, prediction_horizon );
% eph = SimulateSatellite_integerMPC( satelliteFilename, t0_MJD, ...
%     timestep, duration, prediction_horizon );


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

quaternions = eph(:,8:11);
PlotEulerAngles( quaternions );

PlotQuaternionError( simConfig.referenceQuaternion, quaternions)


%%

%PlotMultipleEulerAngles( eph3, eph2, eph )


%%

if simConfig.enableRW
rw_w = eph(:,15:18);
PlotRWMomentum( rw_w )
end

%%


omega = eph(:,12:14);
PlotRotRate( omega )


%%

% qArray = quaternion(eph(:,8:11));
% tp = theaterPlot;
% op = orientationPlotter(tp);
% axis([-1 1 -1 1 -1 1])
% view(3)
% for ephIter = 1:length(eph(:,1))
%     %clf reset
%     plotOrientation(op,qArray(ephIter,:));
%     pause(1e-3)
% end


%%

if simConfig.enablePropulsion
prop_rho = eph(:,25:30);
PlotPropRho( prop_rho )
com = eph(:,31:33);
PlotCoM( com )
end
