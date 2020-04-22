function res = ga_CostFn( U, B_chi, chi_0, H, W, U_int_index )


global propulsionData


U(U_int_index) = propulsionData.maxThrust .* U(U_int_index);


res = (chi_0 + B_chi*U')'*H*(chi_0 + B_chi*U') + U*W*U';



end

