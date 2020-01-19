function [ eph ] = Ephemeris( Y0, t0, stepLength, numSteps, sat_M_mat, sat_I_mat )

stepTimes = zeros( 1, numSteps );

for stepIter = 1 : numSteps
    
   stepTimes( stepIter ) = t0 + stepIter * stepLength;
    
end

opts = odeset('RelTol',1e-12,'AbsTol',1e-14);

[ t, Y ] = ode45( @(t, Y) SatelliteAcceleration(t, Y, sat_M_mat, sat_I_mat), stepTimes, Y0, opts );

eph = [ t, Y ];

end

