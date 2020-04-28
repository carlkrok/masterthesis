function [ F_srp_body, T_srp_body ] = SRP_ForceTorque( sunVec_body, ...
    coeffR, com_struct, ...
    surfaceCenterVectorsNormalVectorsAreas )

global SRP_CONSTANT

p_srp = SRP_CONSTANT / (3*10^8);

F_srp_body = zeros(3,1);

T_srp_body = zeros(3,1);

sunVec_unit = sunVec_body / norm( sunVec_body );

for surface = surfaceCenterVectorsNormalVectorsAreas
    
    surfaceCenter = (surface(1:3)-com_struct);
    surfaceNormal = surface(4:6);
    surfaceArea = surface(7);
    
    dotProd = dot( surfaceNormal, sunVec_unit);
    
    if dotProd > 0
        
        F_surface = -p_srp * coeffR * dotProd * surfaceArea * sunVec_unit;

        F_srp_body = F_srp_body + F_surface;

        T_surface = cross( surfaceCenter, F_surface );

        T_srp_body = T_srp_body + T_surface;
        
    end
    
end


end

