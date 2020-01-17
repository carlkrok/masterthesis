function [ dY ] = SatelliteRotationalAcceleration( t, Y )

v = Y(4:6);

a = TransAcc_Newton( Y(1:3) );

dY = [ v; a ];

end

