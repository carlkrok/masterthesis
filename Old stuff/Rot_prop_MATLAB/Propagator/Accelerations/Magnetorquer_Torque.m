function t_mtq = Magnetorquer_Torque( m, b_earth, t )

t_mtq = cross( m, b_earth);

if t>1e3
    disp('here')
end

end

