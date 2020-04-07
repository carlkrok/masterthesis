function [ refVec_ECEF ] = getPointingVector( sat_ECEF, target_ECEF )

refVec_ECEF = (target_ECEF - sat_ECEF) / norm(target_ECEF - sat_ECEF);

end

