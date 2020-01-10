function [ eph ] = Ephemeris( Y0, stepLength, numSteps )

stepTimes = zeros( 1, numSteps );

for stepIter = 1 : numSteps
    
   stepTimes( stepIter ) = stepIter * stepLength;
    
end


[ t, Y ] = ode45( @SatelliteAcceleration, stepTimes, Y0 );

eph = [ t, Y ];

end

