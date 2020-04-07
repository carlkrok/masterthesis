function R = EulerToRotMat( z, y, x )

R = [ cos(z)*cos(y),    cos(z)*sin(y)*sin(x) - cos(x)*sin(z),   sin(z)*sin(x) + cos(z)*cos(x)*sin(y);
    cos(y)*sin(z),      cos(z)*cos(x) + sin(z)*sin(y)*sin(x),   cos(x)*sin(z)*sin(y) - cos(z)*sin(x);
    -sin(y),            cos(y)*sin(x),                          cos(y)*cos(x)];

R = [ cos(z)*cos(y),    cos(z)*sin(y)*sin(x) - cos(x)*sin(z),   sin(z)*sin(x) + cos(z)*cos(x)*sin(y);
    cos(y)*sin(z),      cos(z)*cos(x) + sin(z)*sin(y)*sin(x),   cos(x)*sin(z)*sin(y) - cos(z)*sin(x);
    -sin(y),            cos(y)*sin(x),                          cos(y)*cos(x)];

end

