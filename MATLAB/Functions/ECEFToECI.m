function [ r_eci, v_eci, a_eci ] = ECEFToECI( utc, r_ecef, v_ecef, a_ecef )


mjd = mjuliandate(utc);

dUT1=deltaUT1(mjd);
pm = polarMotion(mjd)*180/pi; % 'action', 'none'
dCIP = deltaCIP(mjd);

[ r_eci, v_eci, a_eci ] = ecef2eci(utc,r_ecef,v_ecef,a_ecef,'pm',pm,'dCIP',dCIP,'dUT1',dUT1,'dAT',37);



end

