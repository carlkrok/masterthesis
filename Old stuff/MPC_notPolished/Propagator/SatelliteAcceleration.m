function [ dY ] = SatelliteAcceleration( Y, U, t, mjd )

global simConfig
global rwData
global mtqData
global propulsionData
global satData
global missionData
global plotData
global OMEGA_EARTH

%% Extrancting data from state vector

sim_time = (mjd - missionData.mjd0) * 86400;

r_ECI = Y(1:3);
v_ECI = Y(4:6);
q_ECI = quatnormalize(Y(7:10)')';
w_sat_body = Y(11:13);
w_rw_w = Y(14:17);

I_mat_body = [Y(18), Y(21), Y(22); ...
                Y(21), Y(19), Y(23); ...
                Y(22), Y(23), Y(20)];
            
rho_thrusters = Y(24:29);

com_struct = Y(30:32);


%% Extrancting data from input vector

MTQ_Cmd = U(1:3); % zeros(3,1); %

RW_Cmd = U(4:7);

PROP_Cmd = U(8:13);

            
%% Declaring state vector derivative values

disturbance_torques = zeros(3,1);
a_tot = zeros(3,1);

%% Converting to different reference frames

tot_mass = satData.constr.body_mass;
for thrusterIter = 1:length(satData.propulsion.thrusters)
    [ thisCom, thisMass ] = propMassCoM( satData.propulsion.thrusters(thrusterIter).structureDim, rho_thrusters(thrusterIter) );
    tot_mass = tot_mass + thisMass;
end


v_Rel_ECI = v_ECI - cross( [0; 0; OMEGA_EARTH], r_ECI );
r_ECEF = posECIToECEF(t, mjd, r_ECI);
[ latitude, longitude, altitude ] = ECEFToLLA( r_ECEF );

rotMat_ECIToBody = QuaternionToRotMat(q_ECI);
rotMat_BodyToECI = rotMat_ECIToBody';

NED_zAxis = - r_ECI./norm(r_ECI);
NED_yAxis = cross(NED_zAxis,[0;0;1]);
NED_xAxis = cross(NED_yAxis, NED_zAxis);
rotMat_NEDToECI = [ NED_xAxis, NED_yAxis, NED_zAxis ];

r_Body = rotMat_ECIToBody * r_ECI;
v_Rel_Body = rotMat_ECIToBody * v_Rel_ECI;



%% Calculate Controller Values


b_earth_NED = MagneticField( mjd, latitude, longitude, altitude );
b_earth_ECI = rotMat_NEDToECI*b_earth_NED;
b_earth_body = rotMat_ECIToBody * b_earth_ECI;
    
m_mtq_body = MTQ_Cmd;

MTQ_t = Magnetorquer_Torque( m_mtq_body, b_earth_body, t );


wdot_rw_w = RW_Cmd;


F_prop_cmd = PROP_Cmd;

    
[PROP_f_body, PROP_t_body, PROP_rho_dot, PROP_mass_dot ] = ...
    PROP_ForceTorque( com_struct, ...
    satData.propulsion.thrusters, F_prop_cmd, t );

PROP_f_ECI = rotMat_BodyToECI * PROP_f_body;
a_tot = a_tot + PROP_f_ECI ./ tot_mass;
    



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


a_grav = Gravity_Acc( missionData.mu, r_ECI );
a_tot = a_tot + a_grav;

if simConfig.enableJ2
    a_j2 = J2_Acc( r_ECI, missionData.mu, missionData.J2, missionData.Re);
    a_tot = a_tot + a_j2;
end

if simConfig.enableGravityGradient
    t_grav = Gravity_Torque( r_Body, missionData.mu, I_mat_body );
    disturbance_torques = disturbance_torques + t_grav;
    plotData.grav_torque_norm = [ plotData.grav_torque_norm; sim_time, norm(t_grav) ];
end


if simConfig.enableDrag
    [f_drag_body, t_drag_body] = AtmosDrag_ForceTorque( v_Rel_Body, ...
        satData.coeffDrag, com_struct, ...
        satData.surfaceCenterVectorsNormalVectorsAreas );

    f_drag_ECI = rotMat_BodyToECI * f_drag_body;
    a_tot = a_tot + f_drag_ECI ./ tot_mass;

    disturbance_torques = disturbance_torques + t_drag_body;
    plotData.drag_torque_norm = [ plotData.drag_torque_norm; sim_time, norm(t_drag_body) ];

end

if simConfig.enableSRP
    sunVec_ECI = SunVector( mjd );
    if SIGHT_Alg( r_ECI, sunVec_ECI ) 
        sunVec_body = rotMat_ECIToBody * sunVec_ECI;
        [ f_srp_body, t_srp_body ] = SRP_ForceTorque( sunVec_body, ...
            satData.coeffR, com_struct, ...
            satData.surfaceCenterVectorsNormalVectorsAreas );

        f_srp_ECI = rotMat_BodyToECI * f_srp_body;
        a_tot = a_tot + f_srp_ECI ./ tot_mass;

        disturbance_torques = disturbance_torques + t_srp_body;
        plotData.srp_torque_norm = [ plotData.srp_torque_norm; sim_time, norm(t_srp_body) ];

    end
end
    

%% Assembling the state vector derivative 



wdot_sat_body ...
    =  (I_mat_body + ( rwData.A_mat * rwData.I_mat * rwData.A_MPinv_mat ) ) ...
    \ ( ...
        - ( ...
             rwData.A_mat * rwData.I_mat * wdot_rw_w ...
        ) ...
        + disturbance_torques ...
        + MTQ_t ...
        + PROP_t_body ...
        - cross(w_sat_body, ( (I_mat_body + rwData.A_mat * rwData.I_mat * rwData.A_MPinv_mat) * w_sat_body )) ...
        - I_mat_dot_struct * w_sat_body ...
    );


rdot = v_ECI;
vdot = a_tot;
qdot =  T_q(q_ECI)*w_sat_body;

dY = [ rdot; vdot; qdot; wdot_sat_body; wdot_rw_w; ...
    Ixx_dot_struct; Iyy_dot_struct; Izz_dot_struct; ...
    Ixy_dot_struct; Ixz_dot_struct; Iyz_dot_struct; ...
    PROP_rho_dot; com_dot_struct ];

plotData.wdot_sat_body = [ plotData.wdot_sat_body; sim_time, wdot_sat_body' ];
plotData.disturbance_torques_norm = [ plotData.disturbance_torques_norm; sim_time, norm(disturbance_torques) ];

end

