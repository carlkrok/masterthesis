function [ r_eci ] = posECEFToECI( mjd, r_ecef )

dUT1=deltaUT1(mjd);
pm = polarMotion(mjd)*180/pi; % 'action', 'none'
dCIP = deltaCIP(mjd);

mjdconv = mjd+678942;
utc = datevec(mjdconv);

dcm = dcmeci2ecef('IAU-2000/2006',utc,37,dUT1,pm,'dCIP',dCIP);
r_eci = dcm' * r_ecef;


end