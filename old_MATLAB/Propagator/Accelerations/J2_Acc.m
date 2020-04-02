function [ a ] = J2_Acc( rVec, mu, J2, Re )
% Spacecraft Dynamics and Control p.164
zg = [0; 0; 1];
a = ( (3*mu*J2*(Re^2)) / (2*(norm(rVec)^5)) ) *...
    ( ( ( ( (5*(rVec'*zg)^2)/(norm(rVec)^2) ) - 1) * rVec) - 2*(rVec'*zg)*zg );

end

