
run earthParameters;
run initialConditions_Satellite1;

[ r_PQW, v_PQW ] = OrbitalElementsToPQW( MU_EARTH, sat_h, sat_e, sat_T );
R_PQWToECI = RotMat_PQWToECI( sat_i, sat_O, sat_w );
r_ECI = R_PQWToECI * r_PQW;
v_ECI = R_PQWToECI * v_PQW;

Y0 = [ r_ECI; v_ECI; sat_q_ECI'; sat_w_ECI' ];
t0 = 0;

t_sat = SatellitePeriod( MU_EARTH, a );

stepLength = 1;
numSteps = floor( (10*t_sat) / stepLength );

eph = Ephemeris( Y0, t0, stepLength, numSteps );

figure(1)
hold on
axis equal
grid on
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
plot3( eph(:,2), eph(:,3), eph(:,4) )
hold off


qArray = quaternion(eph(:,8:11));
tp = theaterPlot;
op = orientationPlotter(tp);
plotOrientation(op,qArray);
