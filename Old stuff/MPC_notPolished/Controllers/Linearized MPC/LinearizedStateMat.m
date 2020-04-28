function A = LinearizedStateMat( r, mu, q )


A = [ zeros(3,3), eye(3,3), zeros(3,23); ... r_dot
    -( mu / ( norm( r )^3 )) .* eye(3,3), zeros(3,26); ... v_dot
    zeros(4,10), T_q(q), zeros(4,16); ... q_dot
    zeros(3,29); ... omega_dot
    zeros(4,29); ... rw_dot
    zeros(6,29); ... rho_dot
    zeros(6,29)]; ... thrust_dot

end

