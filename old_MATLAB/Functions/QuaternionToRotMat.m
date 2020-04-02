function [ rotMat ] = QuaternionToRotMat( q )

tol = 1e-6;
if abs(norm(q)-1) > tol 
    error('q must be unit quaternion'); 
end

eta = q(1);
S = S_crossProdMat(q(2:4));

rotMat = eye(3,3)  + 2*eta.*S + 2*S^2;

end

