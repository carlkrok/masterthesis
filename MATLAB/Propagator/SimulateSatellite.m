function [ eph ] = SimulateSatellite( satelliteFilename, mjd0, stepTimes, simConfig )

global MU_EARTH
global R_EARTH
global J2_EARTH

run(satelliteFilename);

[r_PQW, v_PQW] = OrbitalElementsToPQW( MU_EARTH, sat.initCond.orb_h, ...
    sat.initCond.orb_e, sat.initCond.orb_T );
R_PQWToECI = RotMat_PQWToECI( sat.initCond.orb_i, sat.initCond.orb_O, ...
    sat.initCond.orb_w );
r_ECI = R_PQWToECI * r_PQW;
v_ECI = R_PQWToECI * v_PQW;

Y0 = [ r_ECI; v_ECI; sat.initCond.q_ECI; sat.initCond.w_body; sat.rw.h ];


if simConfig.enableRW
    if not( sat.rw.exists )
        error("Satellite not configured with Reaction Wheels")
    end
    rwData.A_mat = sat.rw.A_mat;
    rwData.A_MPinv_mat = sat.rw.A_MPinv_mat;
    rwData.I_mat = sat.rw.I_mat;
    rwData.maxTorque = sat.rw.maxTorque;
    rwData.maxMomentum = sat.rw.maxMomentum;
    rwData.maxVelocity = sat.rw.maxVel;
else
    rwData.A_mat = zeros(3,4);
    rwData.A_MPinv_mat = zeros(4,3);
    rwData.I_mat = zeros(3,3);
    rwData.maxTorque = 0;
    rwData.maxMomentum = 0;
    rwData.maxVelocity = 0;
end
    

satData.M_mat = sat.constr.M_mat;
satData.I_mat = sat.constr.I_mat;

missionData.mjd0 = mjd0;
missionData.mu = MU_EARTH;
missionData.Re = R_EARTH;
missionData.J2 = J2_EARTH;


if simConfig.enablePointing
    missionData.pointingTarget_LLA = simConfig.targetLLA; 
    missionData.pointingTarget_ECEF = LLAToECEF( missionData.pointingTarget_LLA );
else
    missionData.pointingTarget_LLA = [0,0,0];
    missionData.pointingTarget_ECEF = [0;0;0];
end

eph = Ephemeris( Y0, stepTimes, missionData, satData, rwData, simConfig );

end

