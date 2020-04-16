function [ dY ] = SatelliteAcceleration( t, Y )

global simConfig
global rwData
global mtqData
global propulsionData
global satData
global missionData

global OMEGA_EARTH

global plotData

mjd = missionData.mjd0 + t/86400;


%% Extrancting data from state vector

r_ECI = Y(1:3);
v_ECI = Y(4:6);
q_ECI = quatnormalize(Y(7:10)')';
w_sat_body = Y(11:13);
w_rw_w = Y(14:17);

I_mat_struct = [Y(18), Y(21), Y(22); ...
                Y(21), Y(19), Y(23); ...
                Y(22), Y(23), Y(20)];
            
rho_thrusters = Y(24:29);

com_struct = Y(30:32);

qError_ECI_integrated = Y(33:36);
            
%% Declaring state vector derivative values

disturbance_torques = zeros(3,1);
a_tot = zeros(3,1);

%% Converting to different reference frames

tot_mass = satData.constr.body_mass;
for thrusterIter = 1:length(satData.propulsion.thrusters)
    [ thisCom, thisMass ] = propMassCoM( satData.propulsion.thrusters(thrusterIter).structureDim, rho_thrusters(thrusterIter) );
    tot_mass = tot_mass + thisMass;
end

[ I_mat_struct_rightEigenvectors, I_mat_body, I_mat_struct_leftEigenvectors] = eig( I_mat_struct );

[eigValuesSorted,indEigValuesSorted] = sort(diag(satData.constr.I_mat_body));
I_mat_body = I_mat_body(indEigValuesSorted,indEigValuesSorted);
I_mat_struct_rightEigenvectors = I_mat_struct_rightEigenvectors(indEigValuesSorted,indEigValuesSorted);
I_mat_struct_leftEigenvectors = I_mat_struct_leftEigenvectors(indEigValuesSorted,indEigValuesSorted);

RotMat_structToBody = I_mat_struct_rightEigenvectors';


v_Rel_ECI = v_ECI - cross( [0; 0; OMEGA_EARTH], r_ECI );
r_ECEF = posECIToECEF(t, missionData.mjd0, r_ECI);
[ latitude, longitude, altitude ] = ECEFToLLA( r_ECEF );

rotMat_ECIToBody = QuaternionToRotMat(q_ECI);
rotMat_BodyToECI = rotMat_ECIToBody';

NED_zAxis = - r_ECI./norm(r_ECI);
NED_yAxis = cross(NED_zAxis,[0;0;1]);
NED_xAxis = cross(NED_yAxis, NED_zAxis);
rotMat_NEDToECI = [ NED_xAxis, NED_yAxis, NED_zAxis ];

r_Body = rotMat_ECIToBody * r_ECI;
v_Rel_Body = rotMat_ECIToBody * v_Rel_ECI;


%% Find error variables

if simConfig.enablePointing
    targetVec_ECEF = getPointingVector(r_ECEF, simConfig.pointingTarget_ECEF);
    targetVec_ECI = posECEFToECI(t, missionData.mjd0, targetVec_ECEF);
    rotMatRef_ECI = getReferenceFrame( targetVec_ECI, v_ECI );
    qRef_ECI = RotMatToQuaternion( rotMatRef_ECI );
else
    qRef_ECI = simConfig.referenceQuaternion;
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


%% Calculate Controller Values

if simConfig.enableMTQ
    b_earth_NED = MagneticField( mjd, latitude, longitude, altitude );
    b_earth_ECI = rotMat_NEDToECI*b_earth_NED;
    %b_earth_body = rotMat_ECIToBody * b_earth_ECI;
    
    plotData.mtq_b = [plotData.mtq_b; [b_earth_ECI;t]'];

    [m_mtq_body, t_pdcmd_sat_mtq, t_pdcan_sat_mtq] = MTQ_PD( qError_ECI, ...
        omegaError_Body, qError_ECI_integrated, b_earth_ECI, I_mat_body, ...
        simConfig.mtq_controller.Kp, simConfig.mtq_controller.Kd, ...
        simConfig.mtq_controller.Ki, mtqData.maxDipoleMoment, t );
    
    plotData.mtq_m = [plotData.mtq_m; [m_mtq_body;t]'];
    plotData.mtq_t_cmd = [plotData.mtq_t_cmd; [t_pdcmd_sat_mtq;t]'];
    plotData.mtq_t_can = [plotData.mtq_t_can; [t_pdcan_sat_mtq;t]'];
    
    %m_mtq_body = MTQ_BDOT( b_earth_body, t, mtqData.maxDipoleMoment );
    
    MTQ_t = Magnetorquer_Torque( m_mtq_body, b_earth_ECI, t );
    
    plotData.mtq_t = [plotData.mtq_t; [MTQ_t;t]'];

    if t>1e3
        disp('here')
    end

else
    MTQ_t = zeros(3,1);
end


if simConfig.enableRW
    wdot_rw_w = RW_NonlinPD( qError_ECI, omegaError_Body, ...
       I_mat_body, w_sat_body, w_rw_w, rwData.A_mat, rwData.A_MPinv_mat, ...
       rwData.I_mat, rwData.maxAcc, rwData.maxVel, MTQ_t, ...
        simConfig.rw_controller.Kp, simConfig.rw_controller.Kd, ...
        RotMat_structToBody );
else 
    wdot_rw_w = zeros(4,1);
end

if simConfig.enablePropulsion
%     PROP_rho_dot = -1e-3*[1, 2, 3, 0.5, 1.5, 2.5]';
%     PROP_mass_dot = 1e-6*PROP_rho_dot;

    [F_prop_cmd, t_ref_sat_body] = PROP_PD( qError_ECI, omegaError_Body, ...
        simConfig.propulsion_controller.Kp, simConfig.propulsion_controller.Kd, ...
        com_struct, RotMat_structToBody, satData.propulsion.thrusters, ...
        satData.propulsion.maxThrust, t);
    
    [PROP_f_body, PROP_t_body, PROP_rho_dot, PROP_mass_dot ] = ...
        PROP_ForceTorque( com_struct, RotMat_structToBody, ...
        satData.propulsion.thrusters, F_prop_cmd, t );
    
    PROP_f_ECI = rotMat_BodyToECI * PROP_f_body;
    a_tot = a_tot + PROP_f_ECI ./ tot_mass;
    
    %PROP_t_body = t_ref_sat_body;
    
else 
    PROP_rho_dot = zeros(6,1);
    PROP_mass_dot = zeros(6,1);
    
    PROP_t_body = zeros(3,1);
end


%% Calculate state vector derivatives

com_dot_struct = zeros(3,1);
for propIerator = 1:length(PROP_mass_dot)
    [ thisPropCom, thisPropMass ] = propMassCoM( ...
        satData.propulsion.thrusters(propIerator).structureDim, ...
        rho_thrusters(propIerator) );
    com_dot_struct = com_dot_struct + (-com_struct * ...
        (satData.constr.body_mass/tot_mass^2) + thisPropCom * (1/tot_mass) ...
        - thisPropCom * (thisPropMass/tot_mass^2) ...
        ) * PROP_mass_dot(propIerator);
end

Ixx_dot_struct = Inertia_i_dot( satData.constr.body_rho, rho_thrusters, ...
    PROP_rho_dot, satData.constr.body_dim, satData.constr.body_boundaries, ...
    com_struct, com_dot_struct, satData.propulsion.thrusters, 1, 2, 3);

Iyy_dot_struct = Inertia_i_dot( satData.constr.body_rho, rho_thrusters, ...
    PROP_rho_dot, satData.constr.body_dim, satData.constr.body_boundaries, ...
    com_struct, com_dot_struct, satData.propulsion.thrusters, 2, 1, 3);

Izz_dot_struct = Inertia_i_dot( satData.constr.body_rho, rho_thrusters, ...
    PROP_rho_dot, satData.constr.body_dim, satData.constr.body_boundaries, ...
    com_struct, com_dot_struct, satData.propulsion.thrusters, 3, 2, 1);

Ixy_dot_struct  = Inertia_ij_dot( satData.constr.body_rho, rho_thrusters, ....
    PROP_rho_dot, satData.constr.body_boundaries, com_struct, com_dot_struct, ...
    satData.propulsion.thrusters, 1, 2, 3 );

Ixz_dot_struct  = Inertia_ij_dot( satData.constr.body_rho, rho_thrusters, ....
    PROP_rho_dot, satData.constr.body_boundaries, com_struct, com_dot_struct, ...
    satData.propulsion.thrusters, 1, 3, 2 );

Iyz_dot_struct  = Inertia_ij_dot( satData.constr.body_rho, rho_thrusters, ....
    PROP_rho_dot, satData.constr.body_boundaries, com_struct, com_dot_struct, ...
    satData.propulsion.thrusters, 2, 3, 1 );

I_mat_dot_struct = [Ixx_dot_struct, Ixy_dot_struct, Ixz_dot_struct; ...
                    Ixy_dot_struct, Iyy_dot_struct, Iyz_dot_struct; ...
                    Ixz_dot_struct, Iyz_dot_struct, Izz_dot_struct];

I_mat_struct_rightEigenvectors_dot = EigenVectorMatrixDot( ...
    I_mat_struct_rightEigenvectors, I_mat_struct_leftEigenvectors, ...
    eigValuesSorted, I_mat_dot_struct );                

RotMat_structToBody_dot = I_mat_struct_rightEigenvectors_dot';

qStruct_ECI = QuaternionSummation( q_ECI, RotMatToQuaternion( RotMat_structToBody_dot' ));
plotData.qStruct_ECI = [plotData.qStruct_ECI, qStruct_ECI];

wdot_body_body = ExtractOmegaFromSkewSymmetric(RotMat_structToBody_dot * RotMat_structToBody');

a_grav = Gravity_Acc( missionData.mu, r_ECI );
a_tot = a_tot + a_grav;

if simConfig.enableJ2
    a_j2 = J2_Acc( r_ECI, missionData.mu, missionData.J2, missionData.Re);
    a_tot = a_tot + a_j2;
end

if simConfig.enableGravityGradient
    t_grav = Gravity_Torque( r_Body, missionData.mu, I_mat_body );
    disturbance_torques = disturbance_torques + t_grav;
end


if simConfig.enableDrag
    [f_drag_body, t_drag_body] = AtmosDrag_ForceTorque( v_Rel_Body, ...
        satData.coeffDrag, com_struct, RotMat_structToBody, ...
        satData.surfaceCenterVectorsNormalVectorsAreas );

    f_drag_ECI = rotMat_BodyToECI * f_drag_body;
    a_tot = a_tot + f_drag_ECI ./ tot_mass;

    disturbance_torques = disturbance_torques + t_drag_body;
end

if simConfig.enableSRP
    sunVec_ECI = SunVector( mjd );
    if SIGHT_Alg( r_ECI, sunVec_ECI ) 
        sunVec_body = rotMat_ECIToBody * sunVec_ECI;
        [ f_srp_body, t_srp_body ] = SRP_ForceTorque( sunVec_body, ...
            satData.coeffR, com_struct, RotMat_structToBody, ...
            satData.surfaceCenterVectorsNormalVectorsAreas );

        f_srp_ECI = rotMat_BodyToECI * f_srp_body;
        a_tot = a_tot + f_srp_ECI ./ tot_mass;

        disturbance_torques = disturbance_torques + t_srp_body;
    end
end
    

%% Assembling the state vector derivative 

% RW_t_body = RotMat_structToBody * rwData.A_mat * rwData.I_mat * ( ...
%         wdot_rw_w ...
%         + rwData.A_MPinv_mat * RotMat_structToBody' * wdot_sat_body ... 
%         - rwData.A_MPinv_mat * RotMat_structToBody' * wdot_body_body );
% 
% wdot_sat_body = I_mat_body \ ( ...
%     - RW_t_body ...
%     + disturbance_torques ...
%     + MTQ_t ...
%     - cross(w_sat_body, ( I_mat_body * w_sat_body )) ...
%     - cross(w_sat_body, ( RotMat_structToBody * (rwData.A_mat * (rwData.I_mat * w_rw_w)) )) ...
%     - cross(wdot_body_body, ( RotMat_structToBody * rwData.A_mat * rwData.I_mat * w_rw_w )) ...
%     + cross(w_sat_body, ( RotMat_structToBody * rwData.A_mat * ...
%         rwData.I_mat * rwData.A_MPinv_mat * RotMat_structToBody' * wdot_body_body)) ...
%     - RotMat_structToBody * I_mat_dot_struct * RotMat_structToBody' * w_sat_body ...
%     );


wdot_sat_body ...
    =  (I_mat_body + (RotMat_structToBody * rwData.A_mat * rwData.I_mat * + rwData.A_MPinv_mat * RotMat_structToBody' ) ) ...
    \ ( ...
        - ( ...
            RotMat_structToBody * rwData.A_mat * rwData.I_mat * ( ...
                wdot_rw_w ...
                - rwData.A_MPinv_mat * RotMat_structToBody' * wdot_body_body ...
            ) ...
        ) ...
        + disturbance_torques ...
        + MTQ_t ...
        + PROP_t_body ...
        - cross(w_sat_body, ( I_mat_body * w_sat_body )) ...
        - cross(w_sat_body, ( RotMat_structToBody * (rwData.A_mat * (rwData.I_mat * w_rw_w)) )) ...
        - cross(wdot_body_body, ( RotMat_structToBody * rwData.A_mat * rwData.I_mat * w_rw_w )) ...
        + cross(w_sat_body, ( RotMat_structToBody * rwData.A_mat * ...
            rwData.I_mat * rwData.A_MPinv_mat * RotMat_structToBody' * wdot_body_body)) ...
        - RotMat_structToBody * I_mat_dot_struct * RotMat_structToBody' * w_sat_body ...
    );


rdot = v_ECI;
vdot = a_tot;
qdot =  T_q(q_ECI)*w_sat_body;

dY = [ rdot; vdot; qdot; wdot_sat_body; wdot_rw_w; ...
    Ixx_dot_struct; Iyy_dot_struct; Izz_dot_struct; ...
    Ixy_dot_struct; Ixz_dot_struct; Iyz_dot_struct; ...
    PROP_rho_dot; com_dot_struct; qError_ECI ];

if any(isnan(dY), 'all')
    warning("NaN found in satellite acceleration vector")
end

if any(norm(wdot_sat_body) > 1000, 'all')
    warning("Too high satellite rotational velocity")
end

end

