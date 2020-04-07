function q = EulerToQuaternion(xAngle,yAngle,zAngle)
% q = EULER2Q(phi,theta,psi) computes the unit quaternions q = [eta eps1 eps2 eps3]
% from Euler angles phi, theta and psi
%
% Author:    Thor I. Fossen
% Date:      8th June 2000
% Revisions: 6 October 2001, T I. Fossen - eta as first element in q 
%            16 April 2019, T. I. Fossen - replaced Sheppard with NASA algorithm
%
% NASA Mission Planning and Analysis Division. "Euler Angles, Quaternions, 
% and Transformation Matrices" (PDF). NASA. Retrieved 12 January 2013.

cy = cos(zAngle * 0.5);
sy = sin(zAngle * 0.5);
cp = cos(yAngle * 0.5);
sp = sin(yAngle * 0.5);
cr = cos(xAngle * 0.5);
sr = sin(xAngle * 0.5);

q = [cy * cp * cr + sy * sp * sr
     cy * cp * sr - sy * sp * cr
     sy * cp * sr + cy * sp * cr
     sy * cp * cr - cy * sp * sr];

q = q/(q'*q);
