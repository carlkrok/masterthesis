global MU_EARTH
global R_EARTH

% All units in SI (rad, meters)

sat_q_ECI = [1, 0, 0, 0];
sat_w_ECI = [0.00001, 0, 0.00001];

sat_T = 0;
sat_i = 0;
sat_O = 0;
sat_w = 0;

sat_r_a = 500 * 10^3 + R_EARTH;
sat_r_p = 500 * 10^3 + R_EARTH;

sat_a = 0.5 * ( sat_r_a + sat_r_p );
sat_e = Eccentricity(sat_r_p, sat_a);
sat_h = MomentumNorm( MU_EARTH, sat_a, sat_e );

