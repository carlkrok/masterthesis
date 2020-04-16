function B = LinearizedMTQ( I_mat, b_earth )
% Linearized satellite omega derivative response to MTQ gain

B = I_mat \ S_crossProdMat(-b_earth);

end

