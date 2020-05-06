function t_mtq = Magnetorquer_Torque( m, b_earth )

t_mtq = S_crossProdMat(-b_earth)*m;


end

