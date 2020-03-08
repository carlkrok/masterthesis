function [ b ] = MagneticField( mjd, latitude, longitude, altitude )

mjdconv = mjd+678942;

[nBx, nBy, nBz] = igrf(mjdconv, latitude, longitude, altitude, 'geodetic');

b = 10^9\[nBx; nBy; nBz];

end

