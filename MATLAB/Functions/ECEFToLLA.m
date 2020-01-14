function [ latitude, longitude, altitude ] = ECEFToLLA( rECEF )

lla = ecef2lla(rECEF', 'WGS84');

latitude = lla(1); % degrees
longitude = lla(2); % degrees
altitude = lla(3); % meters

end

