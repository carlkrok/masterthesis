function [ torque_body, wdot ] = Forces_RW( w, wdot_commanded, rotMat, I_mat, maxTorque, maxMomentum )

numRW = length(w);

torque_body = zeros(3,1);
wdot = zeros(numRW,1);

for rwIter = 1:numRW
   
    if 
    
end

end

