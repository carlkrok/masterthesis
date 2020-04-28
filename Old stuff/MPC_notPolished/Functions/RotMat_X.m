function [ rotMatX] = RotMat_X( alpha )

rotMatX =   [1, 0,          0; ...
            0,  cos(alpha), -sin(alpha); ...
            0,  sin(alpha),  cos(alpha)];

end

