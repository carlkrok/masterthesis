function [ latitude, longitude, altitude ] = ECEFToLLA( rECI )

lla = ecef2lla(rECI, 'WGS84');

latitude = lla(1); % degrees
longitude = lla(2); % degrees
altitude = lla(3); % meters

end

