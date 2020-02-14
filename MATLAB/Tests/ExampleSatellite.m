clear all

satelliteFilename = "Satellite001";

run earthParameters;
run(satelliteFilename);

t0_UTC = [2016 10 31 12 56 37.67];
t0_MJD = mjuliandate(t0_UTC);

sat_period = SatellitePeriod( MU_EARTH, sat.initCond.orb_a );
stepLength = 10; % [s]
numSteps = floor( (1 * sat_period) / stepLength );
stepTimes = zeros( 1, numSteps+1 );
for stepIter = 0 : numSteps
   stepTimes( stepIter+1 ) = (stepIter * stepLength);
end

simConfig.enableJ2 = false;
simConfig.enableRW = false;
simConfig.enableDrag = true;
simConfig.enableSRP = true;

simConfig.enablePointing = false;
simConfig.targetLLA = [63.4184922, 10.4005655, 0]; 

eph = SimulateSatellite( satelliteFilename, t0_MJD, stepTimes, simConfig );


%%


includeEarth = false;
xyz = eph(:,2:4);
PlotOrbit( xyz, includeEarth );


%%


isCircular = true;
xyz = eph(:,2:4);
t = eph(:,1);
PlotOrbitRadiusDeviation( t, xyz, isCircular );

   
%%

quaternions = eph(:,8:11);
PlotEulerAngles( quaternions );


%%


rw_momentum = eph(:,15:18);
PlotRWMomentum( rw_momentum )


%%

% qArray = quaternion(eph(:,8:11));
% tp = theaterPlot;
% op = orientationPlotter(tp);
% plotOrientation(op,qArray);

