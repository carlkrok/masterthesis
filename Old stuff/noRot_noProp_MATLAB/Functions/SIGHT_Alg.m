function sightBool = SIGHT_Alg( r1, r2 )

global R_EARTH

sightBool = false;

tau_min = (norm(r1)^2 - dot(r1, r2)) / ...
    (norm(r1)^2 + norm(r2)^2 - 2*dot(r1, r2));

if tau_min < 0 || tau_min > 1
    sightBool = true;
elseif (1-tau_min)*norm(r1)^2 + dot(r1, r2)*tau_min >= R_EARTH^2
    sightBool = true;
end

end

