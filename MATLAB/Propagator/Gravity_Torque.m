function [ t ] = Gravity_Torque( rBody, mu, I_mat )

t = (3*mu/(norm(rBody)^5)) * ( cross( rBody, (I_mat * rBody) ) );

end

