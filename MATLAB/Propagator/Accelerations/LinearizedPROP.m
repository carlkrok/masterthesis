function B = LinearizedPROP( I_mat, com_struct, thrusters, numThrusters )

B_raw = [ S_crossProdMat( thrusters(1).structArm - com_struct) * ...
    ( - thrusters(1).structureThrustDir )]; 

for thrustIter = 2:numThrusters
    B_raw = [B_raw, S_crossProdMat( thrusters(thrustIter).structArm - com_struct) * ...
    ( - thrusters(thrustIter).structureThrustDir ) ];
end

B = I_mat \ B_raw;

end

