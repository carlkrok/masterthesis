function [ r_ecef, v_ecef, a_ecef ] = ECIToECEF( mjd, r_eci, v_eci, a_eci )

% dcm=dcmeci2ecef(reduction,utc,deltaAT,deltaUT1,polarmotion)
% Direction cosine or transformation matrix

% dcm = dcmeci2ecef('IAU-2000/2006',[2000 1 12 4 52 12.4], )
% January 12, 2000 at 4 hours, 52 minutes, 12.4 seconds 
% deltaAT � Difference between International Atomic Time and UTC: 37  (TAI is ahead of UTC by this amount)
% deltaUT1 � Difference between UTC and Universal Time (UT1): 2019-05-02: 	-0.2 s
% https://www.nist.gov/pml/time-and-frequency-division/atomic-standards/leap-second-and-ut1-utc-information
% polarmotion � Polar displacement of the Earth, in radians, from the motion of the Earth crust, along the x- and y-axes.
% polarmotion=polarMotion(utc)
% [polarmotion,polarmotionError]=polarMotion(utc) 
%name-value 'dCIP' � Adjustment to the location of the Celestial Intermediate Pole (CIP)
% mjd = mjuliandate(2015,12,28)
% dCIP = deltaCIP(mjd)

dUT1=deltaUT1(mjd);
pm = polarMotion(mjd)*180/pi; % 'action', 'none'
dCIP = deltaCIP(mjd);

mjdconv = mjd+678942;
utc = datevec(mjdconv);

[r_ecef,v_ecef,a_ecef] = eci2ecef(utc,r_eci,v_eci,a_eci,'pm',pm,'dCIP',dCIP,'dUT1',dUT1,'dAT',37);



end

