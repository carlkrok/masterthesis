global MU_EARTH
global R_EARTH

% All units in SI (rad, meters, kg, s)



% Satellite dimensions 
sat.constr.dim_x = 0.2263; % Based on HYPSO_ANA_001 
sat.constr.dim_y = 0.100;
sat.constr.dim_z = 0.366;
sat.constr.body_dim = [sat.constr.dim_x; ...
    sat.constr.dim_y; ...
    sat.constr.dim_z];
sat.constr.body_boundaries = [0, sat.constr.dim_x; ...
    0, sat.constr.dim_y; ...
    0, sat.constr.dim_z];


%% Propulsion
% Based on water propulsion from https://steamjet.space

% NanoFEEP https://doi.org/10.1016/j.actaastro.2018.01.012

sat.propulsion.exists = true;
propRho0 = 10^4;
std_isp = 3000;


%sat.propulsion.minThrust = 15 * 10^-3; % N
sat.propulsion.maxThrust = 20 * 10^-6; % 20 * 10^-3; % N
sat.propulsion.power = 0.84;

% Thruster 1 controlling positive X-axis spin
sat.propulsion.thrusters(1).isp = std_isp;
sat.propulsion.thrusters(1).structureDim = ...
    [sat.constr.dim_x/2-0.005, sat.constr.dim_x/2+0.005; ...
    sat.constr.dim_y-0.01, sat.constr.dim_y; ...
    sat.constr.dim_z-0.01, sat.constr.dim_z];
sat.propulsion.thrusters(1).structureThrustDir = [0;1;0];
sat.propulsion.thrusters(1).rho = propRho0;
sat.propulsion.thrusters(1).structArm = [sat.constr.dim_x/2; ...
    sat.constr.dim_y; sat.constr.dim_z];
% Thruster 2 controlling negative X-axis spin
sat.propulsion.thrusters(2).isp = std_isp;
sat.propulsion.thrusters(2).structureDim = ...
    [sat.constr.dim_x/2-0.005, sat.constr.dim_x/2+0.005; ...
    0, 0.01;
    sat.constr.dim_z-0.01, sat.constr.dim_z ];
sat.propulsion.thrusters(2).structureThrustDir = [0;-1;0];
sat.propulsion.thrusters(2).rho = propRho0;
sat.propulsion.thrusters(2).structArm = [sat.constr.dim_x/2; ...
    0; sat.constr.dim_z ];

% Thruster 3 controlling positive Y-axis spin
sat.propulsion.thrusters(3).isp = std_isp;
sat.propulsion.thrusters(3).structureDim = ...
    [0, 0.01; ...
    sat.constr.dim_y/2-0.005, sat.constr.dim_y/2+0.005; ...
    sat.constr.dim_z-0.01, sat.constr.dim_z];
sat.propulsion.thrusters(3).structureThrustDir = [-1;0;0];
sat.propulsion.thrusters(3).rho = propRho0;
sat.propulsion.thrusters(3).structArm = [0; ...
    sat.constr.dim_y/2; sat.constr.dim_z];
% Thruster 4 controlling negative Y-axis spin
sat.propulsion.thrusters(4).isp = std_isp;
sat.propulsion.thrusters(4).structureDim = ...
    [sat.constr.dim_x-0.01, sat.constr.dim_x; ...
    sat.constr.dim_y/2-0.005, sat.constr.dim_y/2+0.005; ...
    sat.constr.dim_z-0.01, sat.constr.dim_z];
sat.propulsion.thrusters(4).structureThrustDir = [1;0;0];
sat.propulsion.thrusters(4).rho = propRho0;
sat.propulsion.thrusters(4).structArm = [sat.constr.dim_x; ...
    sat.constr.dim_y/2; sat.constr.dim_z];

% Thruster 5 controlling positive Z-axis spin
sat.propulsion.thrusters(5).isp = std_isp;
sat.propulsion.thrusters(5).structureDim = ...
    [0, 0.01; ...
    0, 0.01; ...
    sat.constr.dim_z-0.01, sat.constr.dim_z]; % sat.constr.dim_z-0.01, sat.constr.dim_z]; % 
sat.propulsion.thrusters(5).structureThrustDir = [-sqrt(1/2);sqrt(1/2);0];
sat.propulsion.thrusters(5).rho = propRho0;
sat.propulsion.thrusters(5).structArm = [0; ...
    0; sat.constr.dim_z]; % 0; sat.constr.dim_z]; % 
% Thruster 6 controlling negative Z-axis spin
sat.propulsion.thrusters(6).isp = std_isp;
sat.propulsion.thrusters(6).structureDim = ...
    [sat.constr.dim_x-0.01, sat.constr.dim_x; ...
    sat.constr.dim_y-0.01, sat.constr.dim_y; ...
    sat.constr.dim_z-0.01, sat.constr.dim_z]; % sat.constr.dim_z-0.01, sat.constr.dim_z]; % 
sat.propulsion.thrusters(6).structureThrustDir = [-sqrt(1/2);sqrt(1/2);0];
sat.propulsion.thrusters(6).rho = propRho0;
sat.propulsion.thrusters(6).structArm = [sat.constr.dim_x; ...
    sat.constr.dim_y; sat.constr.dim_z]; % sat.constr.dim_y; sat.constr.dim_z]; % 


%% Satellite construction


% Satellite initial mass and inertia data
sat.constr.body_mass = 5.7; % Based on HYPSO_ANA_001 
sat.constr.body_rho = sat.constr.body_mass / (sat.constr.dim_x * sat.constr.dim_y ...
    * sat.constr.dim_z);

[prop_init_com, prop_init_mass] = totPropMassCoM( sat.propulsion.thrusters );
sat.constr.tot_init_mass = sat.constr.body_mass + prop_init_mass;
sat.constr.tot_init_com_struct = [sat.constr.dim_x/2; sat.constr.dim_y/2; sat.constr.dim_z/2];
sat.constr.body_init_com_struct = (1/sat.constr.body_mass) * ...
    (sat.constr.tot_init_com_struct * sat.constr.tot_init_mass - ...
    prop_init_com * prop_init_mass);



sat.constr.Ixx_struct = Inertia_i( sat.constr.body_rho, propRho0*ones(6,1), ...
    sat.constr.body_dim, sat.constr.body_boundaries, ...
    sat.constr.tot_init_com_struct, sat.propulsion.thrusters, ...
    1, 2, 3);

sat.constr.Iyy_struct = Inertia_i( sat.constr.body_rho, propRho0*ones(6,1), ...
    sat.constr.body_dim, sat.constr.body_boundaries, ...
    sat.constr.tot_init_com_struct, sat.propulsion.thrusters, ...
    2, 1, 3);

sat.constr.Izz_struct = Inertia_i( sat.constr.body_rho, propRho0*ones(6,1), ...
    sat.constr.body_dim, sat.constr.body_boundaries, ...
    sat.constr.tot_init_com_struct, sat.propulsion.thrusters, ...
    3, 1, 2);

sat.constr.Ixy_struct = -1 * Inertia_ij( sat.constr.body_rho, propRho0*ones(6,1), ... 
    sat.constr.body_boundaries, sat.constr.tot_init_com_struct, ...
    sat.propulsion.thrusters, 1, 2, 3);
sat.constr.Iyx_struct = sat.constr.Ixy_struct;

sat.constr.Ixz_struct = -1 * Inertia_ij( sat.constr.body_rho, propRho0*ones(6,1), ... 
    sat.constr.body_boundaries, sat.constr.tot_init_com_struct, ...
    sat.propulsion.thrusters, 1, 3, 2);
sat.constr.Izx_struct = sat.constr.Ixz_struct;

sat.constr.Izy_struct = -1 * Inertia_ij( sat.constr.body_rho, propRho0*ones(6,1), ... 
    sat.constr.body_boundaries, sat.constr.tot_init_com_struct, ...
    sat.propulsion.thrusters, 3, 2, 1);
sat.constr.Iyz_struct = sat.constr.Izy_struct;
                    
sat.constr.I_mat_struct = [sat.constr.Ixx_struct, sat.constr.Ixy_struct, sat.constr.Ixz_struct; ...
                            sat.constr.Iyx_struct, sat.constr.Iyy_struct, sat.constr.Iyz_struct; ...
                            sat.constr.Izx_struct, sat.constr.Izy_struct, sat.constr.Izz_struct];

[ sat.constr.RotMat_structToBody, sat.constr.I_mat_body] = eig( sat.constr.I_mat_struct );
[eigValuesSorted,indEigValuesSorted] = sort(diag(sat.constr.I_mat_body));
sat.constr.I_mat_body = sat.constr.I_mat_body(indEigValuesSorted,indEigValuesSorted);
sat.constr.RotMat_structToBody = sat.constr.RotMat_structToBody(indEigValuesSorted,indEigValuesSorted);

sat.constr.coeffDrag = 2.2;
sat.constr.coeffR = 1.0;

sat.constr.surfaceCenterVectorsNormalVectorsAreas = [ ...
    [sat.constr.dim_x;	sat.constr.dim_y/2;	sat.constr.dim_z/2; ...
        1;	0;	0; ...
        sat.constr.dim_y*sat.constr.dim_z], ...
    [0;	sat.constr.dim_y/2;	sat.constr.dim_z/2; ...
        -1;	0;	0; ...
        sat.constr.dim_y*sat.constr.dim_z], ...
    [sat.constr.dim_x/2;	sat.constr.dim_y;	sat.constr.dim_z/2; ...
        0;	1;	0; ...
        sat.constr.dim_x*sat.constr.dim_z], ...
    [sat.constr.dim_x/2;	0;	sat.constr.dim_z/2; ...
        0;	-1;	0; ...
        sat.constr.dim_x*sat.constr.dim_z], ...
    [sat.constr.dim_x/2;	sat.constr.dim_y/2;	sat.constr.dim_z; ...
        0;	0;	1; ...
        sat.constr.dim_x*sat.constr.dim_y], ...
    [sat.constr.dim_x/2;	sat.constr.dim_y/2;	0; ...
        0;	0;	-1; ...
        sat.constr.dim_x*sat.constr.dim_y], ...
    ];


            
%% Satellite initial condition

% Initial orientation quaternion and angular velocity
sat.initCond.q_ECI = EulerToQuaternion(pi/3,pi/8,-pi/5); % [1; 0; 0; 0]; % (0,0,0); % 
sat.initCond.w_body = [0; 0; 0];
% Initial orbit parameters
sat.initCond.orb_T = 45;
sat.initCond.orb_i = 97;
sat.initCond.orb_O = 0;
sat.initCond.orb_w = 0;

% Satellite orbit apoapsis and periapsis
sat.initCond.orb_r_a = 500 * 10^3 + R_EARTH;
sat.initCond.orb_r_p = 500 * 10^3 + R_EARTH;

sat.initCond.orb_a = 0.5 * ( sat.initCond.orb_r_a + sat.initCond.orb_r_p );
sat.initCond.orb_e = Eccentricity( sat.initCond.orb_r_p, sat.initCond.orb_a);
sat.initCond.orb_h = MomentumNorm( MU_EARTH, sat.initCond.orb_a, sat.initCond.orb_e );




%% Magnetorquer

sat.mtq.exists = true;

sat.mtq.maxDipoleMoment = 2*0.42;

sat.mtq.maxPower = 2*0.86;

sat.mtq.powerFactor = 0.86 / 0.42^2;

sat.mtq.resistance = 5^2 / 0.86; % Voltage squared over power

sat.mtq.alpha = 0.42/5; % nS/R = m/V 




%% Reaction Wheels configuration

sat.rw.exists = true;

sat.rw.h = [0; 0; 0; 0]; % Initial momentum values for RWs    

sat.rw.maxVel = 6500 * (pi / 30); % rad / s
sat.rw.maxMomentum = 20*10^-3;
sat.rw.maxTorque = 3.2*10^-3;
sat.rw.mass = 0.137;

sat.rw.maxPower = 3;
sat.rw.idlePower = 0.045;

sat.rw.w = zeros(4,1); % Initial angular velocity values for RWs 

%sat_reaction_wheels based on nanoavionics NA-4RWO-GO-R8, modelled as solid disk 
sat.rw.radius = sqrt( ( sat.rw.maxMomentum ) / ( sat.rw.maxVel * 0.5 * sat.rw.mass ));
sat.rw.I_mat = (sat.rw.mass * sat.rw.radius^2) * 0.5 * eye(4,4);
sat.rw.A_mat =  [sqrt(2/3), -sqrt(2/3),  0,          0;
                sqrt(1/3),  sqrt(1/3),  -sqrt(1/3),  -sqrt(1/3);
                0,          0,          sqrt(2/3),  -sqrt(2/3);]; % Configuration matrix
            
sat.rw.A_MPinv_mat = sat.rw.A_mat'/(sat.rw.A_mat*sat.rw.A_mat');

sat.rw.maxAcc = sat.rw.maxTorque / ((sat.rw.mass * sat.rw.radius^2) * 0.5);
sat.rw.maxVel = sat.rw.maxMomentum / ((sat.rw.mass * sat.rw.radius^2) * 0.5);

sat.rw.efficiency = 1/sat.rw.maxPower * sat.rw.maxVel * sat.rw.maxTorque;

