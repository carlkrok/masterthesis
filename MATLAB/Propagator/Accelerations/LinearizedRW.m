function B = LinearizedRW( I_mat_sys, I_mat_rw, A_mat_rw )

B = -1.* I_mat_sys \ ( A_mat_rw * I_mat_rw );

end

