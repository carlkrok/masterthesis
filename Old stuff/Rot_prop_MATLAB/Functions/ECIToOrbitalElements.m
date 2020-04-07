function [ angularMomentumNorm, inclination, RAAN, eccentricityNorm, perigeeArgument, trueAnomaly ] = ECIToOrbitalElements( mu, rECI , vECI )

% Based on Alg. 4.2 in Curtis, p. 209

rNorm = norm(rECI);

vNorm = norm(vECI);

vRadial = dot(rECI,vECI) / rNorm;

angularMomentum = cross(rECI,vECI);

angularMomentumNorm = norm(angularMomentum);

inclination = arccos(angularMomentum(3)/angularMomentumNorm);

nodeLine = cross([0;0;1],angularMomentum);

nodeLineNorm = norm(nodeLine);

if nodeLine(2) >= 0
    
    RAAN = arccos(nodeLine(1)/nodeLineNorm);
    
else
    
    RAAN = 2*pi - arccos(nodeLine(1)/nodeLineNorm);
    
end

eccentricity = 1/MU_EARTH * ( ( vNorm^2 - mu/rNorm ) * rECI ...
    - rNorm * norm(vRadial) * vECI );

eccentricityNorm = norm(eccentricity);

if eccentricity(3) >= 0
    perigeeArgument = arccos( dot(nodeLine, eccentricity) / ...
        ( nodeLineNorm * eccentricityNorm ) )
    
else
    
    perigeeArgument = 2*pi - arccos( dot(nodeLine, eccentricity) / ...
        ( nodeLineNorm * eccentricityNorm ) )
    
end


if norm(vRadial) >= 0
    
    trueAnomaly = arccos( dot( (eccentricity/eccentricityNorm), ...
        (rECI/rNorm) ) );
    
else
    
    trueAnomaly = 2*pi - arccos( dot( (eccentricity/eccentricityNorm), ...
        (rECI/rNorm) ) );
    
end


end

