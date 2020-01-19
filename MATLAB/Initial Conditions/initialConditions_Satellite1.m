global MU_EARTH
global R_EARTH

% All units in SI (rad, meters, kg, s)

sat_q_ECI = [1, 0, 0, 0];
sat_w_ECI = [0, 0, 0];

sat_T = 0;
sat_i = 75;
sat_O = 0;
sat_w = 0;

sat_r_a = 500 * 10^3 + R_EARTH;
sat_r_p = 500 * 10^3 + R_EARTH;

sat_mass = 10;
sat_M_mat = sat_mass .* ones(3,3);

Ixx = 10;
Iyy = 10;
Izz = 10;
sat_I_mat = [   Ixx,  0,      0; ...
                0,    Iyy,    0; ...
                0,    0,      Izz ];

sat_a = 0.5 * ( sat_r_a + sat_r_p );
sat_e = Eccentricity(sat_r_p, sat_a);
sat_h = MomentumNorm( MU_EARTH, sat_a, sat_e );

