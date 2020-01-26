function [ dY ] = SatelliteAcceleration( t, Y, mu, pointingTarget_ECEF, sat_M_mat, sat_I_mat, rw_A_mat, rw_I_mat, rw_maxTorque, rw_maxMomentum )

r_ECI = Y(1:3);
v_ECI = Y(4:6);
q_ECI = quatnormalize(Y(7:10)')';
sat_w = Y(11:13);
rw_h = Y(14:17);

r_ECEF = posECIToECEF(t, r_ECI);

targetVec_ECEF = getPointingVector(r_ECEF, pointingTarget_ECEF);
targetVec_ECI = posECEFToECI(t, targetVec_ECEF);

refRotMat = [1,0,0;0,0,1;0,-1,0];%getReferenceFrame( targetVec_ECI, v_ECI );
qRef = RotMatToQuaternion( refRotMat );

% Error variables describe rotation from qurrent attitude to reference
% attitude.
qError = QuaternionError(qRef, q_ECI);

if exist('DEBUGMODE') == 1
    unitQuatTol = 1e-6;
    if abs(norm(qRef)-1) > unitQuatTol
        disp("qRef not unit quaternion")
    end
    if abs(norm(q_ECI)-1) > unitQuatTol
        disp("q_ECI not unit quaternion")
    end
    if abs(norm(qError)-1) > unitQuatTol
        disp("qError not unit quaternion")
    end
end

rotMatError = QuaternionToRotMat( qError );
[ xAngleError, yAngleError, zAngleError ] = RotMatToEuler( rotMatError );
w_sat_ref = 0.5 .* eye(3,3) * [ xAngleError; yAngleError; zAngleError ];

a_Keplerian = TransAcc_Newton( mu, Y(1:3) );

rdot = v_ECI;
vdot = a_Keplerian;
qdot =  1*Omega_w(sat_w)*q_ECI;

rw_moment_body = Forces_RW( rw_h, rw_A_mat, rw_I_mat, rw_maxTorque, rw_maxMomentum );

disturbance_torques = zeros(3,1);
control_torques = zeros(3,1);

wdot_sat = (sat_I_mat) \ ( cross(-sat_w, (sat_I_mat * sat_w)) + control_torques + disturbance_torques );
%0.5.*(w_sat_ref - w_sat_ECI);

wdot_rw = zeros(4,1);

dY = [ rdot; vdot; qdot; wdot_sat; wdot_rw ];

end

