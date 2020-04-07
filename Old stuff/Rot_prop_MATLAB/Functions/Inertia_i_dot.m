function I_i_dot  = Inertia_i_dot( rho_body, rho_prop, rho_prop_dot, ...
    sat_dim, sat_boundaries, com, com_dot, thrusters, i_index, j_index, k_index)


I_i_body_dot = rho_body*(sat_dim(i_index)) * ( ...
        sat_dim(k_index) * ( ...
        sat_boundaries(j_index,1)^2 - sat_boundaries(j_index,2)^2 ...
            + 2 * com(j_index) * ( sat_boundaries(j_index,2) - ...
                sat_boundaries(j_index,1) ) ...
        ) ...
        * com_dot(j_index) ...
        + ...
        sat_dim(j_index) * ( ...
        sat_boundaries(k_index,1)^2 - sat_boundaries(k_index,2)^2 ...
            + 2 * com(k_index) * ( sat_boundaries(k_index,2) - ...
                sat_boundaries(k_index,1) ) ...
        ) ...
        * com_dot(k_index) ...
    );

I_i_prop_dot = 0;

for thrusterIter = 1:length(thrusters)

    thisIdot = ( ...
        (1/3) * rho_prop_dot(thrusterIter) * ( ...
            ( ...
                thrusters(thrusterIter).structureDim(i_index,2) ...
                - thrusters(thrusterIter).structureDim(i_index,1) ...
            ) ...
            * ...
            ( ...
                thrusters(thrusterIter).structureDim(k_index,2) ...
                - thrusters(thrusterIter).structureDim(k_index,1) ...
            ) ...
            * ...
            ( ...
                (thrusters(thrusterIter).structureDim(j_index,2) - com(j_index))^3 ...
                - (thrusters(thrusterIter).structureDim(j_index,1) - com(j_index))^3 ...
            ) ...
        ) ...
        + ...
        (1/3) * rho_prop_dot(thrusterIter) * ( ...
            ( ...
                thrusters(thrusterIter).structureDim(i_index,2) ...
                - thrusters(thrusterIter).structureDim(i_index,1) ...
            ) ...
            * ...
            ( ...
                thrusters(thrusterIter).structureDim(j_index,2) ...
                - thrusters(thrusterIter).structureDim(j_index,1) ...
            ) ...
            * ...
            ( ...
                (thrusters(thrusterIter).structureDim(k_index,2) - com(k_index))^3 ...
                - (thrusters(thrusterIter).structureDim(k_index,1) - com(k_index))^3 ...
            ) ...
        ) ...
        + ...
        rho_prop(thrusterIter) * com_dot(j_index) * ( ...
            ( ...
                thrusters(thrusterIter).structureDim(i_index,2) ...
                - thrusters(thrusterIter).structureDim(i_index,1) ...
            ) ...
            * ...
            ( ...
                thrusters(thrusterIter).structureDim(k_index,2) ...
                - thrusters(thrusterIter).structureDim(k_index,1) ...
            ) ...
            * ...
            ( ...
                thrusters(thrusterIter).structureDim(j_index,1)^2 ...
                - thrusters(thrusterIter).structureDim(j_index,2)^2 ...
                + 2 * com(j_index) * ( ...
                thrusters(thrusterIter).structureDim(j_index,2) ...
                - thrusters(thrusterIter).structureDim(j_index,1) ) ...
            ) ...
        ) ...
        + ...
        rho_prop(thrusterIter) * com_dot(k_index) * ( ...
            ( ...
                thrusters(thrusterIter).structureDim(i_index,2) ...
                - thrusters(thrusterIter).structureDim(i_index,1) ...
            ) ...
            * ...
            ( ...
                thrusters(thrusterIter).structureDim(j_index,2) ...
                - thrusters(thrusterIter).structureDim(j_index,1) ...
            ) ...
            * ...
            ( ...
                thrusters(thrusterIter).structureDim(k_index,1)^2 ...
                - thrusters(thrusterIter).structureDim(k_index,2)^2 ...
                + 2 * com(k_index) * ( ...
                thrusters(thrusterIter).structureDim(k_index,2) ...
                - thrusters(thrusterIter).structureDim(k_index,1) ) ...
            ) ...
        ) ...
        );
        
    
    I_i_prop_dot = I_i_prop_dot + thisIdot;

end

I_i_dot = I_i_body_dot + I_i_prop_dot;

end

