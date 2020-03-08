function [ dY ] = SatelliteAcceleration( t, Y, mjd0, mu, Re, J2, ...
    enableJ2, enableDrag, enableSRP, enableGravityGradient, ...
    enablePointing, pointingTarget_ECEF, sat_mass, sat_I_mat, enableRW, ...
    rw_A_mat, rw_A_MPinv_mat, rw_I_mat, rw_maxTorque, ...
    rw_maxMomentum, controller_rw_Kp, controller_rw_Kd, ...
    enableMTQ, mtq_maxDipoleMoment, controller_mtq_Kp, ...
    controller_mtq_Kd, sat_surfaceCenterVectorsAndAreas, ...
    sat_coeffR, sat_coeffDrag )

global OMEGA_EARTH

mjd = mjd0 + t/86400;

r_ECI = Y(1:3);
v_ECI = Y(4:6);
q_ECI = quatnormalize(Y(7:10)')';
w_sat_body = Y(11:13);
rw_h = Y(14:17);

disturbance_torques = zeros(3,1);
control_torques = zeros(3,1);

a_tot = zeros(3,1);

v_Rel_ECI = v_ECI - cross( [0; 0; OMEGA_EARTH], r_ECI );

r_ECEF = posECIToECEF(t, mjd0, r_ECI);

[ latitude, longitude, altitude ] = ECEFToLLA( r_ECEF );

rotMat_ECIToBody = QuaternionToRotMat(q_ECI);
rotMat_BodyToECI = rotMat_ECIToBody';

r_Body = rotMat_ECIToBody * r_ECI;

v_Rel_body = rotMat_ECIToBody * v_Rel_ECI;

if enablePointing
    targetVec_ECEF = getPointingVector(r_ECEF, pointingTarget_ECEF);
    targetVec_ECI = posECEFToECI(t, mjd0, targetVec_ECEF);
    rotMatRef_ECI = getReferenceFrame( targetVec_ECI, v_ECI );
    qRef_ECI = RotMatToQuaternion( rotMatRef_ECI );
else
    qRef_ECI = Y(18:21);
end

% Error variables describe rotation from qurrent attitude to reference
% attitude.
qError_ECI = QuaternionError(qRef_ECI, q_ECI);

omegaRef_Body = zeros(3,1);
omegaError_Body = w_sat_body - omegaRef_Body;

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


a_grav = Gravity_Acc( mu, r_ECI );
a_tot = a_tot + a_grav;

if enableJ2
    a_j2 = J2_Acc( r_ECI, mu, J2, Re);
    a_tot = a_tot + a_j2;
end

if enableGravityGradient
    t_grav = Gravity_Torque( r_Body, mu, sat_I_mat );
    disturbance_torques = disturbance_torques + t_grav;
end

if enableMTQ
    b_earth_ECI = MagneticField( mjd, latitude, longitude, altitude );
    b_earth_body = rotMat_ECIToBody * b_earth_ECI;

    m_mtq_body = MTQ_PD( qError_ECI, omegaError_Body, b_earth_body, sat_I_mat, ...
        controller_mtq_Kp, controller_mtq_Kd, mtq_maxDipoleMoment, ...
        w_sat_body, rw_A_mat, rw_h );
    
    %m_mtq_body = MTQ_BDOT( b_earth_body, t, mtq_maxDipoleMoment );
    
    t_mtq_body = Magnetorquer_Torque( m_mtq_body, b_earth_body );
    
    control_torques = control_torques + t_mtq_body;

else
    t_MTQ = zeros(3,1);
end


if enableRW
    rw_commandedTorque = RW_NonlinPD( qError_ECI, omegaError_Body, ...
        sat_I_mat, w_sat_body, rw_h, rw_A_mat, rw_A_MPinv_mat, ...
        rw_maxTorque, rw_maxMomentum, t_MTQ, ...
        controller_rw_Kp, controller_rw_Kd );
    
     [rw_hdot, t_rw] = RW_Torque( w_sat_body, rw_commandedTorque, rw_h, ...
         rw_A_mat, rw_A_MPinv_mat );
    
    control_torques = control_torques + t_rw;
    
else 
    rw_hdot = zeros(4,1);
end

if enableDrag
    [f_drag_body, t_drag_body] = AtmosDrag_ForceTorque( v_Rel_body, sat_coeffDrag, sat_surfaceCenterVectorsAndAreas );

    f_drag_ECI = rotMat_BodyToECI * f_drag_body;
    a_tot = a_tot + f_drag_ECI ./ sat_mass;

    disturbance_torques = disturbance_torques + t_drag_body;
end

if enableSRP
    sunVec_ECI = SunVector( mjd );
    if SIGHT_Alg( r_ECI, sunVec_ECI ) 
        sunVec_body = rotMat_ECIToBody * sunVec_ECI;
        [ f_srp_body, t_srp_body ] = SRP_ForceTorque( sunVec_body, sat_coeffR, sat_surfaceCenterVectorsAndAreas );

        f_srp_ECI = rotMat_BodyToECI * f_srp_body;
        a_tot = a_tot + f_srp_ECI ./ sat_mass;

        disturbance_torques = disturbance_torques + t_srp_body;
    end
end
    

wdot_sat = (sat_I_mat) \ ( ...
    -1 .* (S_crossProdMat(w_sat_body) * sat_I_mat * w_sat_body) ...
    + disturbance_torques ...
    + control_torques ...
    );

rdot = v_ECI;
vdot = a_tot;
qdot =  T_q(q_ECI)*w_sat_body;

dY = [ rdot; vdot; qdot; wdot_sat; rw_hdot; zeros(4,1) ];

if any(isnan(dY), 'all')
    warning("NaN found in satellite acceleration vector")
end

if any(norm(wdot_sat) > 1000, 'all')
    warning("Too high satellite rotational velocity")
end

end

