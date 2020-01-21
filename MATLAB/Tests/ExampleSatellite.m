
run earthParameters;
run Satellite001;

t0_UTC = [2016 10 31 12 56 37.67];
t0_MJD = mjuliandate(t0_UTC);

[ r_PQW, v_PQW ] = OrbitalElementsToPQW( MU_EARTH, sat_h, sat_e, sat_T );
R_PQWToECI = RotMat_PQWToECI( sat_i, sat_O, sat_w );
r_ECI = R_PQWToECI * r_PQW;
v_ECI = R_PQWToECI * v_PQW;

sat_period = SatellitePeriod( MU_EARTH, sat_a );

stepLength = 0.001; % [s]
numSteps = floor( (1 * sat_period) / stepLength );
stepTimes = zeros( 1, numSteps );

qRefTraj = ones( 4, numSteps );
quarterNum = floor(numSteps/4);
qRefTraj(:, 1:quarterNum) = sat_q_ECI .* ones( 4, quarterNum );
qRefTraj(:, quarterNum+1:2*quarterNum) = EulerToQuaternion(pi/2,0,0) .* ones( 4, quarterNum );
qRefTraj(:, 2*quarterNum+1:3*quarterNum) = sat_q_ECI .* ones( 4, quarterNum );
qRefTraj(:, 3*quarterNum+1:end) = EulerToQuaternion(0,pi/2,0) .* ones( 4, numSteps - 3*quarterNum );

Y0 = [ r_ECI; v_ECI; sat_q_ECI; sat_w_ECI; sat_rw_vel ];
t0 = t0_MJD;

for stepIter = 1 : numSteps
   stepTimes( stepIter ) = t0 + stepIter * stepLength;
end

rwData.rotMat = sat_rw_rotMat;
rwData.maxVel = sat_rw_maxVel;
rwData.maxAcc = sat_rw_maxAcc;
rwData.inertia = sat_rw_inertia;

satData.M_mat = sat_M_mat;
satData.I_mat = sat_I_mat;

eph = Ephemeris( Y0, t0, stepTimes, qRefTraj, MU_EARTH, satData, rwData );

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
