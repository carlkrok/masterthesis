function [ eph ] = Ephemeris( Y0, t0, stepLength, numSteps )

stepTimes = zeros( 1, numSteps );

for stepIter = 1 : numSteps
    
   stepTimes( stepIter ) = t0 + stepIter * stepLength;
    
end


[ t, Y ] = ode45( @SatelliteTranslationalAcceleration, stepTimes, Y0 );

eph = [ t, Y ];

end

