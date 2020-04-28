function [com, mass] = totPropMassCoM( thrusters )

mass = 0;
com = [0;0;0];

for thrusterIter = 1:length(thrusters)
    
    [ thisCom, thisMass ] = propMassCoM( thrusters(thrusterIter).structureDim, thrusters(thrusterIter).rho );
   
    mass = mass + thisMass;
    com = com + thisCom.*thisMass;
end

com = com./mass;

end

