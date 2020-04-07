function [ F_drag_body, T_drag_body ] = AtmosDrag_ForceTorque( vRel_body, coeffDrag, surfaceCenterVectorsAndAreas )

vRel_unitVec = vRel_body / norm(vRel_body);

rho = getAtmosphericDensity( 500 * 10^3 );

F_drag_body = zeros(3,1);

T_drag_body = zeros(3,1);

for surface = surfaceCenterVectorsAndAreas
    
    surfaceCenter = surface(1:3);
    surfaceNormal = surfaceCenter/norm(surfaceCenter);
    surfaceArea = surface(4);
    
    dotProd = dot( surfaceNormal, vRel_unitVec);
    
    if dotProd > 0
        
        F_surface = (0.5 * coeffDrag * rho * norm(vRel_body)^2 * dotProd * surfaceArea * vRel_unitVec);

        F_drag_body = F_drag_body + F_surface;

        T_surface = cross( surfaceCenter, F_surface );

        T_drag_body = T_drag_body + T_surface;
        
    end
    
end

end 

