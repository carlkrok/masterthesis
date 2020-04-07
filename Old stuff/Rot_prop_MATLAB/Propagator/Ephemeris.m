function [ eph ] = Ephemeris( Y0, stepTimes )

opts = odeset('RelTol',1e-5,'AbsTol',1e-6,'OutputFcn',@odeprog,'Events',@odeabort);
[t, Y] = ode45( @(t, Y) SatelliteAcceleration(t, Y), stepTimes, Y0, opts );

unitQuatTol = 1e-6;
for stepIter = 1:length(Y(:,1))
    thisNorm = norm(Y(stepIter,7:10));
    if abs(thisNorm - 1) > unitQuatTol
        %fprintf('Non-unit quaternion found in step %i of ephemeris, length %d - will be normalized.\n',stepIter,thisNorm )
        Y(stepIter,7:10) = quatnormalize( Y(stepIter,7:10) );
    end
end

eph = [ t, Y ];

end

