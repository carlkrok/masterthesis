
run earthParameters;
run Satellite001;

t0_UTC = [2016 10 31 12 56 37.67];
t0_MJD = mjuliandate(t0_UTC);

[ r_PQW, v_PQW ] = OrbitalElementsToPQW( MU_EARTH, sat_h, sat_e, sat_T );
R_PQWToECI = RotMat_PQWToECI( sat_i, sat_O, sat_w );
r_ECI = R_PQWToECI * r_PQW;
v_ECI = R_PQWToECI * v_PQW;

Y0 = [ r_ECI; v_ECI; sat_q_ECI'; sat_w_ECI' ];
t0 = t0_MJD;

sat_period = SatellitePeriod( MU_EARTH, sat_a );

stepLength = 1; % [s]
numSteps = floor( (3 * sat_period) / stepLength );

eph = Ephemeris( Y0, t0, stepLength, numSteps, sat_M_mat, sat_I_mat );

%%

figure(1)
hold on
axis equal
grid on
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
plot3( eph(:,2), eph(:,3), eph(:,4) )
hold off



radiusDeviations = zeros(numSteps, 1);
for stepIter = 1:numSteps
    radiusDeviations(stepIter) = norm( [eph(stepIter,2), eph(stepIter,3), eph(stepIter,4)] );
    if sat_r_a == sat_r_p
        radiusDeviations(stepIter) = radiusDeviations(stepIter) - sat_r_a;
    end 
end
figure(2)
hold on
grid on
xlabel('Time [s]')
ylabel('Deviation [m]')
plot(eph(:,1), radiusDeviations)

   
%%

qArray = quaternion(eph(:,8:11));
tp = theaterPlot;
op = orientationPlotter(tp);
plotOrientation(op,qArray);
