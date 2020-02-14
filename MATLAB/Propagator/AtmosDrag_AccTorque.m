function [ a, t ] = AtmosDrag_AccTorque( satMass, satVel, satPos, coeffDrag, area, omergaEarth )

vRel = satVel - cross( [0; 0; omergaEarth], satPos );

a = 0.5 * (1 / satMass) * coeffDrag * area * rho * norm(vRel)^2 * (vRel / norm(vRel));

end 

