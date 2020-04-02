function [ rotMatZ] = RotMat_Z( alpha )

rotMatZ =   [cos(alpha),    -sin(alpha),    0; ...
            sin(alpha),     cos(alpha),     0; ...
            0,              0,              1];

end

