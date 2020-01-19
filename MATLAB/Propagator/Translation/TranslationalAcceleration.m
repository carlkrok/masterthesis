function [ a ] = TranslationalAcceleration( mu, sat_M, t, Y )

a = TransAcc_Newton( mu, Y(1:3) );

end

