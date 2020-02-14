function [ dY ] = SatelliteAcceleration( t, Y, mjd0, mu, Re, J2, enableJ2, enableDrag, enableSRP, enablePointing, pointingTarget_ECEF, sat_M_mat, sat_I_mat, enableRW, rw_A_mat, rw_A_MPinv_mat, rw_I_mat, rw_maxTorque, rw_maxMomentum )

mjd = mjd0 + t/86400;

r_ECI = Y(1:3);
v_ECI = Y(4:6);
q_ECI = quatnormalize(Y(7:10)')';
w_sat_body = Y(11:13);
rw_h = Y(14:17);

disturbance_torques = zeros(3,1);
control_torques = zeros(3,1);

a_tot = zeros(3,1);

r_ECEF = posECIToECEF(t, mjd0, r_ECI);

[ latitude, longitude, altitude ] = ECEFToLLA( r_ECEF );

rotMat_ECIToBody = QuaternionToRotMat(q_ECI);
r_Body = rotMat_ECIToBody * r_ECI;

if enablePointing
    targetVec_ECEF = getPointingVector(r_ECEF, pointingTarget_ECEF);
    targetVec_ECI = posECEFToECI(t, mjd0, targetVec_ECEF);
    rotMatRef_ECI = getReferenceFrame( targetVec_ECI, v_ECI );
    qRef_ECI = RotMatToQuaternion( rotMatRef_ECI );
else
    rotMatRef_ECI = [1,0,0;0,0,1;0,-1,0];
    qRef_ECI = RotMatToQuaternion( rotMatRef_ECI );
end

% Error variables describe rotation from qurrent attitude to reference
% attitude.
qError_ECI = QuaternionError(qRef_ECI, q_ECI);


unitQuatTol = 1e-6;
if abs(norm(qRef_ECI)-1) > unitQuatTol
    disp("qRef not unit quaternion")
end
if abs(norm(q_ECI)-1) > unitQuatTol
    disp("q_ECI not unit quaternion")
end
if abs(norm(qError_ECI)-1) > unitQuatTol
    disp("qError not unit quaternion")
end


% rotMatError_ECI = QuaternionToRotMat( qError_ECI );
% [ xAngleError, yAngleError, zAngleError ] = RotMatToEuler( rotMatError_ECI );
% w_sat_ref = 0.5 .* eye(3,3) * [ xAngleError; yAngleError; zAngleError ];

a_grav = Gravity_Acc( mu, r_ECI );
a_tot = a_tot + a_grav;

if enableJ2
    a_j2 = J2_Acc( r_ECI, mu, J2, Re);
    a_tot = a_tot + a_j2;
end

t_grav = Gravity_Torque( r_Body, mu, sat_I_mat );
disturbance_torques = disturbance_torques + t_grav;

b_earth = MagneticField( mjd, latitude, longitude, altitude );

if enableRW
    rw_appliedTorque = zeros(3,1);
    [t_rw, rw_hdot] = RW_Torque( w_sat_body, rw_appliedTorque, rw_h, rw_A_mat, rw_A_MPinv_mat, rw_maxTorque, rw_maxMomentum );
    control_torques = control_torques + t_rw;
else 
    rw_hdot = zeros(4,1);
end

wdot_sat = (sat_I_mat) \ ( -1 .* cross(w_sat_body, (sat_I_mat * w_sat_body)) + control_torques + disturbance_torques );

rdot = v_ECI;
vdot = a_tot;
qdot =  T_q(q_ECI)*w_sat_body;

dY = [ rdot; vdot; qdot; wdot_sat; rw_hdot ];

end

