function [ dY ] = SatelliteTotalAcceleration( t, Y )

v = Y(4:6);

a = SatelliteTranslationalAcceleration( t, Y );

q = Y(7:10);
w = Y(11:13);
qdot = T_q( q ) * w;
wdot = zeros(3,1);

dY = [ v; a; qdot; wdot ];

end

