global MU_EARTH
global R_EARTH

% All units in SI (rad, meters, kg, s)

% Initial orientation quaternion and angular velocity
sat_q_ECI = [1; 0; 0; 0];
sat_w_body = [0; pi/12; 0];

% Initial orbit parameters
sat_T = 0;
sat_i = 75;
sat_O = 0;
sat_w = 0;

sat_r_a = 500 * 10^3 + R_EARTH;
sat_r_p = 500 * 10^3 + R_EARTH;

sat_a = 0.5 * ( sat_r_a + sat_r_p );
sat_e = Eccentricity(sat_r_p, sat_a);
sat_h = MomentumNorm( MU_EARTH, sat_a, sat_e );

% Satellite dimensions 
sat_dim_x = 0.2263; % Based on HYPSO_ANA_001 
sat_dim_y = 0.100;
sat_dim_z = 0.366;

% Satellite initial mass and inertia data
sat_mass = 5.7; % Based on HYPSO_ANA_001 
sat_M_mat = sat_mass .* eye(3,3);
sat_cg_from_co = [0; 0; 0];

% Assume solid uniform cuboid
Ixx = (1/12) * sat_mass * ( sat_dim_y^2 + sat_dim_z^2 );
Iyy = (1/12) * sat_mass * ( sat_dim_x^2 + sat_dim_z^2 );
Izz = (1/12) * sat_mass * ( sat_dim_x^2 + sat_dim_y^2 );
sat_I_mat = [   Ixx,  0,      0; ...
                0,    Iyy,    0; ...
                0,    0,      Izz ];


sat_rw_maxVel = 6500 * (pi / 30); % rad / s
sat_rw_maxMomentum = 20*10^-3;
sat_rw_maxTorque = 3.2*10^-3;
sat_rw_mass = 0.137;
sat_rw_w = zeros(4,1); % Initial angular velocity values for RWs     
sat_rw_h = [0; 0.5*sat_rw_maxMomentum; 0; 0.5*sat_rw_maxMomentum]; % Initial momentum values for RWs  
%sat_reaction_wheels based on nanoavionics NA-4RWO-GO-R8, modelled as solid disk 
sat_rw_radius = sqrt( ( sat_rw_maxMomentum ) / ( sat_rw_maxVel * 0.5 * sat_rw_mass ));
sat_rw_I_mat = (sat_rw_mass * sat_rw_radius^2) .* [1/4,     0,      0;
                                                    0,      1/4,    0;
                                                    0,      0,      1/2];
sat_rw_A_mat =  [sqrt(2/3), -sqrt(2/3),  0,          0;
                sqrt(1/3),  sqrt(1/3),  -sqrt(1/3),  -sqrt(1/3);
                0,          0,          sqrt(2/3),  -sqrt(2/3);]; % Configuration matrix
sat_rw_A_MPinv_mat = sat_rw_A_mat'/(sat_rw_A_mat*sat_rw_A_mat');

%sat_rw_maxAcc = sat_rw_maxTorque / sat_rw_I_mat(3,3);

% Rotation matrices rotate from RW frame to body frame
%  sat_rw_num = 4;
%  sat_rw_alpha = 0.397904493026690; % Angle of RW's relative to XY plane
%  sat_rw_rotMat = zeros(3,3,sat_rw_num);
%  sat_rw_rotMat(:,:,1) = RotMat_X( (pi/2) - sat_rw_alpha );
%  sat_rw_rotMat(:,:,2) = RotMat_X( sat_rw_alpha - (pi/2) );
%  sat_rw_rotMat(:,:,3) = RotMat_Y( (pi/2) - sat_rw_alpha );
%  sat_rw_rotMat(:,:,4) = RotMat_Y( sat_rw_alpha - (pi/2) );
