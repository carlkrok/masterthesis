function [ a ] = Acceleration_Newton( r )

global MU_EARTH
a = - MU_EARTH * r / ( norm( r )^3 );

end

