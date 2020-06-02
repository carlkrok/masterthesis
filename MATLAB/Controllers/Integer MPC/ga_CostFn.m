function res = ga_CostFn( U, c, D, e, U_int_index )


global propulsionData

U(U_int_index) = propulsionData.thrustRange(U(U_int_index));


res = 2 * c * U' + U * D * U' + e * abs(U');



end

