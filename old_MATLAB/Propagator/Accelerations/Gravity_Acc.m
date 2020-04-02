function [ a ] = Gravity_Acc( mu, r )

a = - ( mu / ( norm( r )^3 )) .* r;

end

