function [ eph ] = SimulateSatellite( satelliteFilename, mjd0, ...
    stepTimes )

global MU_EARTH
global R_EARTH
global J2_EARTH

global simConfig
global rwData
global mtqData
global propulsionData
global satData
global missionData

run(satelliteFilename);

[r_PQW, v_PQW] = OrbitalElementsToPQW( MU_EARTH, sat.initCond.orb_h, ...
    sat.initCond.orb_e, sat.initCond.orb_T );
R_PQWToECI = RotMat_PQWToECI( sat.initCond.orb_i, sat.initCond.orb_O, ...
    sat.initCond.orb_w );
r_ECI = R_PQWToECI * r_PQW;
v_ECI = R_PQWToECI * v_PQW;

Y0 = [ r_ECI; v_ECI; sat.initCond.q_ECI; sat.initCond.w_body; sat.rw.w; ...
    sat.constr.Ixx_struct; sat.constr.Iyy_struct; sat.constr.Izz_struct; ...
    sat.constr.Ixy_struct; sat.constr.Ixz_struct; sat.constr.Iyz_struct; ...
    sat.propulsion.thrusters(1).rho; ...
    sat.propulsion.thrusters(2).rho; ...
    sat.propulsion.thrusters(3).rho; ...
    sat.propulsion.thrusters(4).rho; ...
    sat.propulsion.thrusters(5).rho; ...
    sat.propulsion.thrusters(6).rho; ...
    sat.constr.body_init_com_struct;
    zeros(4,1)];


if not( sat.rw.exists )
    error("Satellite not configured with Reaction Wheels")
end
rwData.A_mat = sat.rw.A_mat;
rwData.A_MPinv_mat = sat.rw.A_MPinv_mat;
rwData.I_mat = sat.rw.I_mat;
rwData.maxAcc = sat.rw.maxAcc;
rwData.maxVel = sat.rw.maxVel;
rwData.maxVelocity = sat.rw.maxVel;

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
else
    propulsionData.minThrust = 0;
    propulsionData.maxThrust = 0;
end

satData.constr = sat.constr;
satData.propulsion = sat.propulsion;
satData.surfaceCenterVectorsNormalVectorsAreas = sat.constr.surfaceCenterVectorsNormalVectorsAreas;
satData.coeffR = sat.constr.coeffR;
satData.coeffDrag = sat.constr.coeffDrag;

missionData.mjd0 = mjd0;
missionData.mu = MU_EARTH;
missionData.Re = R_EARTH;
missionData.J2 = J2_EARTH;

eph = Ephemeris( Y0, stepTimes );

end

