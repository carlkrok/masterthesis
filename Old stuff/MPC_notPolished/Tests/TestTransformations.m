
% asbCubeSatBlockLib
% open('asbCubeSatSimulationProject.sltx')

r0_ECI = [8000*10^3; 0; 0];
v0_ECI = [0; 8000; 0];
a0_ECI = [-1000; -1000; 0];

UTC = [2016 10 31 12 56 37.67];
MJD = mjuliandate(UTC);

%%

[r_ECEF, v_ECEF, a_ECEF] = ECIToECEF( MJD, r0_ECI, v0_ECI, a0_ECI );

[ r_ECI, v_ECI, a_ECI ] = ECEFToECI( MJD, r_ECEF, v_ECEF, a_ECEF );

rDiff = r_ECI - r0_ECI

vDiff = v_ECI - v0_ECI

aDiff = a_ECI - a0_ECI

%%

[ latitude, longitude, altitude ] = ECEFToLLA( r_ECEF );
[ r2_ECEF ] = LLAToECEF( [latitude, longitude, altitude] );
r2_ECEF = r2_ECEF';

rDiff_1 = r2_ECEF - r_ECEF

[ latitude, longitude, altitude ] = ECEFToLLA( r2_ECEF );
[ r3_ECEF ] = LLAToECEF( [latitude, longitude, altitude] );
r3_ECEF = r3_ECEF';

rDiff_2 = r3_ECEF - r2_ECEF

[ latitude, longitude, altitude ] = ECEFToLLA( r3_ECEF );
[ r4_ECEF ] = LLAToECEF( [latitude, longitude, altitude] );
r4_ECEF = r4_ECEF';

rDiff_3 = r4_ECEF - r3_ECEF

%%


[r_ECEF_fn1, v_ECEF_fn1, a_ECEF_fn1] = ECIToECEF( MJD, r0_ECI, v0_ECI, a0_ECI );
r_ECEF_fn2 = posECIToECEF( MJD, r0_ECI );
r_ECI_fn2 = posECEFToECI( MJD, r_ECEF_fn2 );

fnDiff = r_ECEF_fn2 - r_ECEF_fn1
fn2Acc = r_ECI_fn2 - r0_ECI



