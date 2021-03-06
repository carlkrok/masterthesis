function [ com, mass ] = propMassCoM( boundaries, rho )


thisV = ...
    (boundaries(1,2)- ...
    boundaries(1,1)) ...
    * ...
    (boundaries(2,2)- ...
    boundaries(2,1)) ...
    * ...
    (boundaries(3,2)- ...
    boundaries(3,1));

mass = thisV*rho;


com = [(boundaries(1,1) + boundaries(1,2))/2; ...
    (boundaries(2,1) + boundaries(2,2))/2; ...
    (boundaries(3,1) + boundaries(3,2))/2];

end

