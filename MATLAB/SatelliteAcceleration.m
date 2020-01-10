function [ dY ] = SatelliteAcceleration( t, Y )

v = Y(4:6);

a = Acceleration_Newton( Y(1:3) );

dY = [ v; a ];

end

