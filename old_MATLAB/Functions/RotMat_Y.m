function [ rotMatY] = RotMat_Y( alpha )

rotMatY =   [cos(alpha),    0,  sin(alpha); ...
            0,              1,  0; ...
            -sin(alpha),    0,  cos(alpha)];

end

