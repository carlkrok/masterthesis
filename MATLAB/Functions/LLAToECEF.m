function [ rECEF ] = LLAToECEF( LatLonAlt )

rECEF = lla2ecef( LatLonAlt, 'WGS84')';

end

