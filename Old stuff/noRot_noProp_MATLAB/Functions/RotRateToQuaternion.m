function [ q ] = RotRateToQuaternion( omega )
    
    if norm(omega) == 0
        q = [1;0;0;0];
    else
        q = [ cos( norm(omega) / 2 );
        (omega / norm(omega)) * sin( norm(omega) / 2 )];
    end

end

