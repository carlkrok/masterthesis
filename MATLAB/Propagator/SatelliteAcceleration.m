function [ dY ] = SatelliteAcceleration( t, Y, mu, pointingTarget_ECEF, sat_M_mat, sat_I_mat, sat_rw_rotMat )

r0_ECI = Y(1:3);
v0_ECI = Y(4:6);
q0_ECI = Y(7:10);
w0_sat_ECI = Y(11:13);
w0_rw = Y(14:17);

a_Keplerian = TransAcc_Newton( mu, Y(1:3) );

r0_ECEF = posECIToECEF(t, r0_ECI);

targetVec_ECEF = getPointingVector(r0_ECEF, pointingTarget_ECEF);
targetVec_ECI = posECEFToECI(t, targetVec_ECEF);

refRotMat = getReferenceFrame( targetVec_ECI, v0_ECI );
qRef = rotm2quat( refRotMat );

qError = QuaternionError(qRef, q0_ECI);



rdot = v0_ECI;
vdot = a_Keplerian;
qdot =  qError; %0.5*Omega_w(w0_ECI)*q0_ECI;
wdot_sat = w0_sat_ECI;
wdot_rw = zeros(4,1);

dY = [ rdot; vdot; qdot; wdot_sat; wdot_rw ];

end

