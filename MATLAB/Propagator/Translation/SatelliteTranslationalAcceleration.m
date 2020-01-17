function [ a ] = SatelliteTranslationalAcceleration( t, Y )

a = TransAcc_Newton( Y(1:3) );

end

