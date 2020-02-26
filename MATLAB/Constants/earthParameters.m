global R_EARTH
global MU_EARTH
global J2_EARTH
global OMEGA_EARTH
global DENSITY_MAX_500KM
global DENSITY_MIN_500KM
global SRP_CONSTANT
global AU_TO_KM
global JD_MJD_DIFF

R_EARTH = 6.378137 * 10^6; % m
MU_EARTH = 3.986004418 * 10^(14); % m^3 (solar s)^-2
J2_EARTH = -0.0010826267; %
OMEGA_EARTH = 0.0000729211585530; % rad / (solar a)

DENSITY_MAX_500KM = 2.042 * 10^(-12); % kg / m^3, from Vallado, Harris-Priester model
DENSITY_MIN_500KM = 3.916 * 10^(-13); % kg / m^3, from Vallado, Harris-Priester model

SRP_CONSTANT = 1367; % W/M^2, from Vallado (Baker) 

AU_TO_KM = 149597870; % km

JD_MJD_DIFF = 2400000.5;
