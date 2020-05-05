function A = LinearizedStateMat( q )


A = [ zeros(4,4), T_q(q), zeros(4,10); ... q_dot
    zeros(3,17); ... omega_dot
    zeros(4,17); ... rw_dot
    zeros(6,17)]; ... rho_dot

end

