function [ dY ] = SatelliteAcceleration( t, Y, sat_M_mat, sat_I_mat )

global MU_EARTH

v = Y(4:6);

a = TransAcc_Newton( MU_EARTH, Y(1:3) );

q = Y(7:10);
w = Y(11:13);
qdot = T_q( q ) * w;

wdot = zeros(3,1);

dY = [ v; a; qdot; wdot ];

end

