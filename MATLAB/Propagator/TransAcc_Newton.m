function [ a ] = TransAcc_Newton( mu, r )

a = - ( mu / ( norm( r )^3 )) .* r;

end

