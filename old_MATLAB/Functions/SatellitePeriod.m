function [ t ] = SatellitePeriod( mu, a )

t = 2 * ( pi / sqrt( mu ) ) * a^(3/2);

end

