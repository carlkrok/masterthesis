function [ eph ] = Ephemeris( Y0, stepTimes, missionData, satData, rwData  )

opts = odeset('RelTol',1e-6,'AbsTol',1e-6);

[ t, Y ] = ode45( @(t, Y) SatelliteAcceleration(t, Y, missionData.mu, missionData.pointingTarget_ECEF, satData.M_mat, satData.I_mat, rwData.A_mat, rwData.I_mat, rwData.maxTorque, rwData.maxMomentum ), stepTimes, Y0, opts );

for stepIter = 1:length(stepTimes)
    if norm(Y(stepIter,7:10)) > 1e-6
        Y(stepIter,7:10) = quatnormalize( Y(stepIter,7:10) );
    end
end

eph = [ t, Y ];

end

