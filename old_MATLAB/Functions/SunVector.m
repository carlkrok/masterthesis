function sunVec = SunVector( MJD )

global AU_TO_KM
global JD_MJD_DIFF;

% Algorithm 29 in Vallado

JD = MJD + JD_MJD_DIFF;

T = (JD - 2451545.0)/36525;

lambda_M = 280.460 + 36000.771 * T; % [deg]

M = 357.5277233 + 35999.05034 * T; % [deg]

lambda_ecliptic = lambda_M + 1.914666471 * sind( M ) ...
    + 0.019994643 * sind( 2*M );

r = 1.000140612 - 0.016708617 * cosd( M ) - 0.000139589 * cosd( 2*M );

epsilon = 23.439291 - 0.0130042 * T;

sunVecAU =    [r * cosd( lambda_ecliptic );
            r * cosd( epsilon ) * sind( lambda_ecliptic );
            r * sind( epsilon ) * sind( lambda_ecliptic )];

sunVec = sunVecAU.*AU_TO_KM;
        
end

