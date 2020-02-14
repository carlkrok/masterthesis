function [ b ] = MagneticField( mjd, latitude, longitude, altitude )

mjdconv = mjd+678942;

[Bx, By, Bz] = igrf(mjdconv, latitude, longitude, altitude, 'geodetic');

b = [Bx; By; Bz];

end

