function [ eph ] = SimulateSatellite( satelliteFilename, mjd0, ...
    stepTimes, simConfig )

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

Y0 = [ r_ECI; v_ECI; sat.initCond.q_ECI; sat.initCond.w_body; sat.rw.h; ...
    simConfig.referenceQuaternion ];


if simConfig.enableRW
    if not( sat.rw.exists )
        error("Satellite not configured with Reaction Wheels")
    end
    rwData.A_mat = sat.rw.A_mat;
    rwData.A_MPinv_mat = sat.rw.A_MPinv_mat;
    rwData.I_mat = sat.rw.I_mat;
    rwData.maxAcc = sat.rw.maxAcc;
    rwData.maxVel = sat.rw.maxVel;
    rwData.maxVelocity = sat.rw.maxVel;
else
    rwData.A_mat = zeros(3,4);
    rwData.A_MPinv_mat = zeros(4,3);
    rwData.I_mat = zeros(3,3);
    rwData.maxAcc = 0;
    rwData.maxVel = 0;
    rwData.maxVelocity = 0;
end

if simConfig.enableMTQ
    if not( sat.mtq.exists )
        error("Satellite not configured with Magnetorquer")
    end
    mtqData.maxDipoleMoment = sat.mtq.maxDipoleMoment;
else
    mtqData.maxDipoleMoment = 0;
end

if simConfig.enablePropulsion
    if not( sat.propulsion.exists )
        error("Satellite not configured with Propulsion system")
    end
    propulsionData.minThrust = sat.propulsion.minThrust;
    propulsionData.maxThrust = sat.propulsion.maxThrust;
    propulsionData.xArm = sat.propulsion.xArm;
    propulsionData.yArm = sat.propulsion.yArm;
    propulsionData.zArm = sat.propulsion.zArm;
else
    propulsionData.minThrust = 0;
    propulsionData.maxThrust = 0;
    propulsionData.xArm = 0;
    propulsionData.yArm = 0;
    propulsionData.zArm = 0;
end

satData.M_mat = sat.constr.M_mat;
satData.mass = sat.constr.mass;
satData.I_mat = sat.constr.I_mat;
satData.surfaceCenterVectorsAndAreas = sat.constr.surfaceCenterVectorsAndAreas;
satData.coeffR = sat.constr.coeffR;
satData.coeffDrag = sat.constr.coeffDrag;

missionData.mjd0 = mjd0;
missionData.mu = MU_EARTH;
missionData.Re = R_EARTH;
missionData.J2 = J2_EARTH;

eph = Ephemeris( Y0, stepTimes, missionData, satData, rwData, mtqData, ...
    propulsionData, simConfig );

end

