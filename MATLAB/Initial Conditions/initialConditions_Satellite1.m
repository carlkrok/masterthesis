global MU_EARTH
global R_EARTH

% All units in SI (rad, meters)

T = 0;
i = 0;
O = 0;
w = 0;

r_a = 500 * 10^3 + R_EARTH;
r_p = 500 * 10^3 + R_EARTH;

a = 0.5 * ( r_a + r_p );
e = Eccentricity(r_p, a);
h = MomentumNorm( MU_EARTH, a, e );

