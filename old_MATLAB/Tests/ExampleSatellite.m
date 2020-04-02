clear all

satelliteFilename = "Satellite001";

run earthParameters;
run(satelliteFilename);

global mtq_m
global mtq_t
global mtq_t_cmd
global mtq_b
mtq_b = zeros(4,1)';
mtq_t_cmd = zeros(4,1)';
mtq_m = zeros(4,1)';
mtq_t = zeros(4,1)';

t0_UTC = [2005 10 31 12 56 37.67];
t0_MJD = mjuliandate(t0_UTC);

sat_period = SatellitePeriod( MU_EARTH, sat.initCond.orb_a );
stepLength = 1; % [s]
numSteps = floor( (0.05 * sat_period) / stepLength );
stepTimes = zeros( 1, numSteps+1 );
for stepIter = 0 : numSteps
   stepTimes( stepIter+1 ) = (stepIter * stepLength);
end

simConfig.enableJ2 = false;
simConfig.enableDrag = false;
simConfig.enableSRP = false;
simConfig.enableGravityGradient = false;

simConfig.enableRW = true;
%simConfig.enableRW = false;
simConfig.rw_controller.Kp = 0.18; % Hz
simConfig.rw_controller.Kd = 0.6;

%simConfig.enableMTQ = true;
simConfig.enableMTQ = false;
simConfig.mtq_controller.Kp = 1e-4; 
simConfig.mtq_controller.Kd = 1e-2;

simConfig.enablePropulsion = false;
simConfig.propulsion_controller.Kp = 0.1; % Hz
simConfig.propulsion_controller.Kd = 0.01;

simConfig.enablePointing = false;
simConfig.pointingTarget_LLA = [63.4184922, 10.4005655, 0];
simConfig.pointingTarget_ECEF = LLAToECEF( simConfig.pointingTarget_LLA );
simConfig.referenceQuaternion = [1; 0; 0; 0];

eph = SimulateSatellite( satelliteFilename, t0_MJD, stepTimes, ...
    simConfig );


%%


% includeEarth = false;
% xyz = eph(:,2:4);
% PlotOrbit( xyz, includeEarth );


%%


% isCircular = true;
% xyz = eph(:,2:4);
% t = eph(:,1);
% PlotOrbitRadiusDeviation( t, xyz, isCircular );

   
%%

quaternions = eph(:,8:11);
PlotEulerAngles( quaternions );


qRef = [1; 0; 0; 0];
PlotQuaternionError( qRef, quaternions)


%%

%PlotMultipleEulerAngles( eph3, eph2, eph )


%%

if simConfig.enableRW
rw_momentum = eph(:,15:18);
PlotRWMomentum( rw_momentum )
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
%     pause(.01)
% end

%%

%startNum = 500;
if simConfig.enableMTQ
    PlotMTQDipole( mtq_m )
    %PlotMTQDipole( mtq_b )
    %PlotMTQTorque( mtq_t )
    PlotTwoTorques( mtq_t, mtq_t_cmd )
    %Plot3DMTQTorquePD( mtq_t(startNum:end,:), mtq_t_cmd(startNum:end,:).*1e-1,...
    %    mtq_m(startNum:end,:).*1e-6, mtq_b(startNum:end,:).*1e-2 )
end






