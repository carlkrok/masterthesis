function [ dY ] = SatelliteAcceleration( t, Y, mu, qRefTraj, sat_M_mat, sat_I_mat, sat_rw_rotMat )

v = Y(4:6);

a = TransAcc_Newton( mu, Y(1:3) );

q = Y(7:10);
w = Y(11:13);
qdot = 0.5*Omega_w(w)*q;

[wdot_sat, wdot_rw] = AngAcc_RW( t, Y, qRefTraj, sat_I_mat, sat_rw_rotMat );

dY = [ v; a; qdot; wdot_sat; wdot_rw ];

end

