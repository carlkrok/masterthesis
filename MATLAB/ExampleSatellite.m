
run earthParameters;
run initialConditions_Satellite1;

[ r_PQW, v_PQW ] = OrbitalElementsToPQW( MU_EARTH, h, e, T );
R_PQWToECI = RotMat_PQWToECI( i, O, w );
r_ECI = R_PQWToECI * r_PQW;
v_ECI = R_PQWToECI * v_PQW;

Y0 = [ r_ECI; v_ECI ];
t0 = 0;

t_sat = SatellitePeriod( MU_EARTH, a );

stepLength = 100;
numSteps = floor( t_sat / stepLength );

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
