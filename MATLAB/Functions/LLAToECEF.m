function [ rECEF ] = LLAToECEF( latitude, longitude, altitude )

rECEF = lla2ecef( [latitude, longitude, altitude], 'WGS84');

end

