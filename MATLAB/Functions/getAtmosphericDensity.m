function dens = getAtmosphericDensity( r )

if norm(r) > 520 || norm(r) < 480
    error("Satellite outside of Atmospheric Density range")
end

global DENSITY_MAX_500KM
global DENSITY_MIN_500KM

dens = 0.5*(DENSITY_MAX_500KM + DENSITY_MIN_500KM);

end

