global MU_EARTH
global R_EARTH

% All units in SI (rad, meters, kg, s)

%% Satellite construction

% Satellite dimensions 
sat.constr.dim_x = 0.2263; % Based on HYPSO_ANA_001 
sat.constr.dim_y = 0.100;
sat.constr.dim_z = 0.366;

% Satellite initial mass and inertia data
sat.constr.mass = 5.7; % Based on HYPSO_ANA_001 
sat.constr.M_mat = sat.constr.mass .* eye(3,3);
sat.constr.cg_from_co = [0; 0; 0];

% Assume solid uniform cuboid
sat.constr.Ixx = (1/12) * sat.constr.mass * ( sat.constr.dim_y^2 + sat.constr.dim_z^2 );
sat.constr.Iyy = (1/12) * sat.constr.mass * ( sat.constr.dim_x^2 + sat.constr.dim_z^2 );
sat.constr.Izz = (1/12) * sat.constr.mass * ( sat.constr.dim_x^2 + sat.constr.dim_y^2 );
sat.constr.I_mat = [   sat.constr.Ixx,    0,                  0; ...
                        0,                  sat.constr.Iyy,    0; ...
                        0,                  0,                  sat.constr.Izz ];
                    
sat.constr.coeffDrag = 2.2;
sat.constr.coeffR = 1.0;

sat.constr.surfaceCenterVectorsAndAreas = [ ...
    [sat.constr.dim_x;  0;                  0; ...
    sat.constr.dim_y*sat.constr.dim_z], ...
    [-sat.constr.dim_x; 0;                  0; ...
    sat.constr.dim_y*sat.constr.dim_z], ...
    [0;                 sat.constr.dim_y;   0; ...
    sat.constr.dim_x*sat.constr.dim_z], ...
    [0;                 -sat.constr.dim_y;  0; ...
    sat.constr.dim_x*sat.constr.dim_z], ...
    [0;                 0;                  sat.constr.dim_z; ...
    sat.constr.dim_x*sat.constr.dim_y], ...
    [0;                 0;                  -sat.constr.dim_z; ...
    sat.constr.dim_x*sat.constr.dim_y]];


            
%% Satellite initial condition

% Initial orientation quaternion and angular velocity
sat.initCond.q_ECI = EulerToQuaternion(pi/2.1,pi/2.1,pi/2.1); % (0,0,0); %
sat.initCond.w_body = [0; 0; 0];
% Initial orbit parameters
sat.initCond.orb_T = 0;
sat.initCond.orb_i = 97;
sat.initCond.orb_O = 0;
sat.initCond.orb_w = 0;

% Satellite orbit apoapsis and periapsis
sat.initCond.orb_r_a = 500 * 10^3 + R_EARTH;
sat.initCond.orb_r_p = 500 * 10^3 + R_EARTH;

sat.initCond.orb_a = 0.5 * ( sat.initCond.orb_r_a + sat.initCond.orb_r_p );
sat.initCond.orb_e = Eccentricity( sat.initCond.orb_r_p, sat.initCond.orb_a);
sat.initCond.orb_h = MomentumNorm( MU_EARTH, sat.initCond.orb_a, sat.initCond.orb_e );


%% Propulsion
% Based on water propulsion from https://steamjet.space

sat.propulsion.exists = true;

sat.propulsion.minThrust = 15 * 10^3; % N
sat.propulsion.maxThrust = 20 * 10^3; % N

sat.propulsion.xArm = sat.constr.dim_x/2;
sat.propulsion.yArm = sat.constr.dim_y/2;
sat.propulsion.zArm = sat.constr.dim_z/2;


%% Magnetorquer

sat.mtq.exists = true;

sat.mtq.maxDipoleMoment = 2*0.42;

sat.mtq.resistance = 5^2 / 0.86; % Voltage squared over power

sat.mtq.alpha = 0.42/5; % nS/R = m/V 




%% Reaction Wheels configuration

sat.rw.exists = true;

sat.rw.h = [0; 0; 0; 0]; % Initial momentum values for RWs    

sat.rw.maxVel = 6500 * (pi / 30); % rad / s
sat.rw.maxMomentum = 20*10^-3;
sat.rw.maxTorque = 3.2*10^-3;
sat.rw.mass = 0.137;
sat.rw.w = zeros(4,1); % Initial angular velocity values for RWs 

%sat_reaction_wheels based on nanoavionics NA-4RWO-GO-R8, modelled as solid disk 
sat.rw.radius = sqrt( ( sat.rw.maxMomentum ) / ( sat.rw.maxVel * 0.5 * sat.rw.mass ));
sat.rw.I_mat = (sat.rw.mass * sat.rw.radius^2) .* [1/4,     0,      0;
                                                    0,      1/4,    0;
                                                    0,      0,      1/2];
sat.rw.A_mat =  [sqrt(2/3), -sqrt(2/3),  0,          0;
                sqrt(1/3),  sqrt(1/3),  -sqrt(1/3),  -sqrt(1/3);
                0,          0,          sqrt(2/3),  -sqrt(2/3);]; % Configuration matrix
sat.rw.A_MPinv_mat = sat.rw.A_mat'/(sat.rw.A_mat*sat.rw.A_mat');

sat1.rw.maxAcc = sat1.rw.maxTorque / sat.rw.mass;
sat1.rw.maxVel = sat1.rw.maxMomentum / sat.rw.mass;

% Rotation matrices rotate from RW frame to body frame
%  sat1.rw.num = 4;
%  sat1.rw.alpha = 0.397904493026690; % Angle of RW's relative to XY plane
%  sat1.rw.rotMat = zeros(3,3,sat1.rw.num);
%  sat1.rw.rotMat(:,:,1) = RotMat_X( (pi/2) - sat1.rw.alpha );
%  sat1.rw.rotMat(:,:,2) = RotMat_X( sat1.rw.alpha - (pi/2) );
%  sat1.rw.rotMat(:,:,3) = RotMat_Y( (pi/2) - sat1.rw.alpha );
%  sat1.rw.rotMat(:,:,4) = RotMat_Y( sat1.rw.alpha - (pi/2) );
