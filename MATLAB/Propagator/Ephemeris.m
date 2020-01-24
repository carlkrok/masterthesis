function [ eph ] = Ephemeris( Y0, stepTimes, missionData, satData, rwData  )

opts = odeset('RelTol',1e-12,'AbsTol',1e-14);

[ t, Y ] = ode45( @(t, Y) SatelliteAcceleration(t, Y, missionData.mu, missionData.pointingTarget_ECEF, satData.M_mat, satData.I_mat, rwData.rotMat ), stepTimes, Y0, opts );

eph = [ t, Y ];

end

