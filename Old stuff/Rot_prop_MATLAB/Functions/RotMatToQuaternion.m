function [ q ] = RotMatToQuaternion( rotMat )

q = rotm2quat( rotMat )';

end

