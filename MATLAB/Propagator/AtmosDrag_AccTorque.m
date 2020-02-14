function [ F_drag, T_drag ] = AtmosDrag_ForceTorque( satVel, satPos, coeffDrag, area, omegaEarth )

rho = getAtmosphericDensity( );

vRel = satVel - cross( [0; 0; omegaEarth], satPos );

F_drag = 0.5 * coeffDrag * area * rho * norm(vRel)^2 * (vRel / norm(vRel));

center_drag = [0;0;1];

T_drag = cross( center_drag, F_drag );

end 

