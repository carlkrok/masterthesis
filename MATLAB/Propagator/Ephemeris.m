function [ eph ] = Ephemeris( Y0, stepTimes, missionData, satData, rwData, simConfig )

%opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
%t, Y = ode45( @(t, Y) SatelliteAcceleration(t, Y, missionData.mjd0, ...
%    missionData.mu, missionData.Re, missionData.J2, ...
%    missionData.pointingTarget_ECEF, satData.M_mat, ...
%    satData.I_mat, rwData.A_mat, rwData.A_MPinv_mat, rwData.I_mat, ...
%    rwData.maxTorque, rwData.maxMomentum ), stepTimes, Y0, opts );

[t, Y] = RK4( Y0, stepTimes, missionData, satData, rwData, simConfig );

unitQuatTol = 1e-6;
for stepIter = 1:length(stepTimes)
    thisNorm = norm(Y(stepIter,7:10));
    if abs(thisNorm - 1) > unitQuatTol
        fprintf('Non-unit quaternion found in step %i of ephemeris, length %d - will be normalized.\n',stepIter,thisNorm )
        Y(stepIter,7:10) = quatnormalize( Y(stepIter,7:10) );
    end
end

eph = [ t, Y ];

end

