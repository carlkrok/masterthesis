function [ hNorm ] = MomentumNorm( mu, a, e )

hNorm = sqrt( mu * a * ( 1 - e^2));

end

