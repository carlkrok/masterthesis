function [ F_srp, T_srp ] = SRP_ForceTorque( )

F_srp = [1;0;0];

center_srp = [0;0;1];

T_srp = cross( center_srp, F_srp );


end

