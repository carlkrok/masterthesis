function [ eph ] = Ephemeris( Y0, t0, stepTimes, qRefTraj, mu, satData, rwData  )

opts = odeset('RelTol',1e-12,'AbsTol',1e-14);

[ t, Y ] = ode45( @(t, Y) SatelliteAcceleration(t, Y, mu, qRefTraj, satData.M_mat, satData.I_mat, rwData.rotMat ), stepTimes, Y0, opts );

eph = [ t, Y ];

end

