function [ b ] = MagneticField( mjd, latitude, longitude, altitude )

mjdconv = mjd+678942;

% Altitude is converted from m to km
[nBx, nBy, nBz] = igrf(mjdconv, latitude, longitude, altitude/10^3, 'geodetic');

b = 10^9.\[nBx; nBy; nBz];

end

