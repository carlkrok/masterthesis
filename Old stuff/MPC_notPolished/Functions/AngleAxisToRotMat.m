function R = AngleAxisToRotMat( angle, axis )

R = cos(angle) * eye(3,3) + sin(angle) * S_crossProdMat(axis) + ...
    (1 - cos(angle)) * (axis * axis');

end

