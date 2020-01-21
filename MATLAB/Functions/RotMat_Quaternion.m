function [ rotMat ] = RotMat_Quaternion( q )

rotMat = eye(3,3)  + 2*q(1).*S_crossProdMat(q(2:4)) + 2*(S_crossProdMat(q(2:4)))^2;

end

