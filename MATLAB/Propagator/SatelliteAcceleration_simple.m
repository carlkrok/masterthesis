function [ dY ] = SatelliteAcceleration_simple( Y, U, ...
    I_mat_body, com_struct, tot_mass, b_earth_eci )

global rwData
global satData
global missionData

t = 0;

%% Extrancting data from state vector

r_ECI = Y(1:3);
v_ECI = Y(4:6);
q_ECI = quatnormalize(Y(7:10)')';
w_sat_body = Y(11:13);
w_rw_w = Y(14:17);
rho_thrusters = Y(18:23);


%% Extrancting data from input vector

MTQ_Cmd = U(1:3);

RW_Cmd = U(4:7);

PROP_Cmd = U(8:13);

            
%% Declaring state vector derivative values

a_tot = zeros(3,1);

%% Converting to different reference frames


rotMat_ECIToBody = QuaternionToRotMat(q_ECI);
rotMat_BodyToECI = rotMat_ECIToBody';

%% Calculate Controller Values

b_earth_body = rotMat_ECIToBody * b_earth_eci;

MTQ_t = Magnetorquer_Torque( MTQ_Cmd, b_earth_body, t );

wdot_rw_w = RW_Cmd;
    
[PROP_f_body, PROP_t_body, PROP_rho_dot, PROP_mass_dot ] = ...
    PROP_ForceTorque( com_struct, ...
    satData.propulsion.thrusters, PROP_Cmd, t );

PROP_f_ECI = rotMat_BodyToECI * PROP_f_body;
%a_tot = a_tot + PROP_f_ECI ./ tot_mass;
    

%% Calculate state vector derivatives


a_grav = Gravity_Acc( missionData.mu, r_ECI );
a_tot = a_tot + a_grav;


%% Assembling the state vector derivative 



wdot_sat_body ...
    =  (I_mat_body + ( rwData.A_mat * rwData.I_mat * rwData.A_MPinv_mat ) ) ...
    * ( ...
        - ( ...
             rwData.A_mat * rwData.I_mat * wdot_rw_w ...
        ) ... 
        + MTQ_t ... 
        + PROP_t_body ...
        - cross(w_sat_body, ( (I_mat_body + rwData.A_mat * rwData.I_mat * rwData.A_MPinv_mat) * w_sat_body )) ...
    );


rdot = v_ECI;
vdot = a_tot;
qdot =  T_q(q_ECI)*w_sat_body;

dY = [ rdot; vdot; qdot; wdot_sat_body; wdot_rw_w; ...
    PROP_rho_dot ];


end

