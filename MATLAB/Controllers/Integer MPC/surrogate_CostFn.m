function res = surrogate_CostFn( U, c, D, e, U_int_index, A_u,b_u )


global propulsionData

U(U_int_index) = propulsionData.maxThrust .* U(U_int_index);


res.Fval = 2 * c * U' + U * D * U' + e * abs(U');

res.Ineq = A_u * U' - b_u;



end
