function A = LinearizedStateMat( q, numPropellant )


A = [ zeros(4,4), T_q(q), zeros(4,4), zeros(4,numPropellant); ... q_dot
    zeros(3,11+numPropellant); ... omega_dot
    zeros(4,11+numPropellant); ... rw_dot
    zeros(numPropellant,11+numPropellant)]; ... rho_dot


end

