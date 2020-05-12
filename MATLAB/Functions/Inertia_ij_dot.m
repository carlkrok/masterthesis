function I_ij_dot  = Inertia_ij_dot( rho_body, rho_prop, rho_prop_dot, ...
    sat_boundaries, com, com_dot, thrusters, i_index, j_index, ...
    k_index, numPropellant )

global simConfig

I_i_body_dot = (rho_body/2) * ( ...
        ( ...
            (sat_boundaries(i_index,1) - sat_boundaries(i_index,2)) ...
            * ...
            (sat_boundaries(k_index,2) - sat_boundaries(k_index,1)) ...
            * ...
            ( ...
                sat_boundaries(j_index,2)^2 - sat_boundaries(j_index,1)^2 ...
                + 2 * ( sat_boundaries(j_index,1) - ...
                    sat_boundaries(j_index,2) ) ...
            ) ...
            * com_dot(i_index) ...
        ) ...
        + ...
        ( ...
            (sat_boundaries(j_index,1) - sat_boundaries(j_index,2)) ...
            * ...
            (sat_boundaries(k_index,2) - sat_boundaries(k_index,1)) ...
            * ...
            ( ...
                sat_boundaries(i_index,2)^2 - sat_boundaries(i_index,1)^2 ...
                + 2 * ( sat_boundaries(i_index,1) - ...
                    sat_boundaries(i_index,2) ) ...
            ) ...
            * com_dot(i_index) ...
        ) ...
        ); 

I_i_prop_dot = 0;

if simConfig.enablePropulsionInertia
    for thrusterIter = 1:numPropellant

        thisIdot = (1/4) * rho_prop_dot(thrusterIter) * ( ...
            ( ...
                thrusters(thrusterIter).structureDim(i_index,2)^2 ...
                - thrusters(thrusterIter).structureDim(i_index,1)^2 ...
                + 2 * ( thrusters(thrusterIter).structureDim(i_index,1) - ...
                thrusters(thrusterIter).structureDim(i_index,2) ) * com(i_index) ...
            ) ...
            * ...
            ( ...
                thrusters(thrusterIter).structureDim(j_index,2)^2 ...
                - thrusters(thrusterIter).structureDim(j_index,1)^2 ...
                + 2 * ( thrusters(thrusterIter).structureDim(j_index,1) - ...
                thrusters(thrusterIter).structureDim(j_index,2) ) * com(j_index) ...
            ) ...
            * ( thrusters(thrusterIter).structureDim(k_index,2) ...
            - thrusters(thrusterIter).structureDim(k_index,1) ) ...
            ) ...
            + ...
            (rho_prop(thrusterIter)/2) * ...
                ( ...
                    ( thrusters(thrusterIter).structureDim(i_index,1) ...
                    - thrusters(thrusterIter).structureDim(i_index,2) ) ...
                * ...
                    ( thrusters(thrusterIter).structureDim(k_index,2) ...
                    - thrusters(thrusterIter).structureDim(k_index,1) ) ...
                * ...
                ( ...
                    thrusters(thrusterIter).structureDim(j_index,2)^2 ...
                    - thrusters(thrusterIter).structureDim(j_index,1)^2 ...
                    + 2 * ( thrusters(thrusterIter).structureDim(j_index,1) - ...
                    thrusters(thrusterIter).structureDim(j_index,2) ) * com(j_index) ...
                ) ...
                ) * com_dot(i_index) ...
                + ...
                ( ...
                    ( thrusters(thrusterIter).structureDim(j_index,1) ...
                    - thrusters(thrusterIter).structureDim(j_index,2) ) ...
                * ...
                    ( thrusters(thrusterIter).structureDim(k_index,2) ...
                    - thrusters(thrusterIter).structureDim(k_index,1) ) ...
                * ...
                ( ...
                    thrusters(thrusterIter).structureDim(i_index,2)^2 ...
                    - thrusters(thrusterIter).structureDim(i_index,1)^2 ...
                    + 2 * ( thrusters(thrusterIter).structureDim(i_index,1) - ...
                    thrusters(thrusterIter).structureDim(i_index,2) ) * com(i_index) ...
                ) ...
                ) * com_dot(j_index); 

        I_i_prop_dot = I_i_prop_dot + thisIdot;

    end
end

I_ij_dot = I_i_body_dot + I_i_prop_dot;

end

