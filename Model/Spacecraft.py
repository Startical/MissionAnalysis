import numpy as np
from CommonTools.Timeseries import Timeseries
from CommonTools.Constants import EARTH_MU, EARTH_RADIUS, EARTH_W
from OrbitTools.FrameTransformations import FrameTransformations as frames
from astropy.coordinates import CartesianRepresentation, EarthLocation, ITRS, GCRS
from astropy.time import Time
import astropy.units as u

import CommonTools.Mathlib as Math
import CommonTools.Datelib as Dates




def orbitalPeriod(sma):
    """Compute the orbital period"""

    T = 2*np.pi*np.sqrt(sma**3/EARTH_MU)

    return T


def compute_lof_sc_guidance_with_constant_bias(xyz_j2000_ts, q_lof_sc = [0,0,0,1]):

    q_j2000_lof_ts = Timeseries('q_j2000_lof',xyz_j2000_ts.refTime)
    q_j2000_sc_ts = Timeseries('q_j2000_sc',xyz_j2000_ts.refTime)

    for i in range(len(xyz_j2000_ts.time)):

        xyz = xyz_j2000_ts.data[i]

        if i == 0:
            xyz_dot = [1,0,0]

            if np.linalg.norm(np.cross(xyz,xyz_dot))<1e-3:
                xyz_dot = [0,1,0]

        else:
            xyz_prev = xyz_j2000_ts.data[i-1]
            dt = xyz_j2000_ts.time[i] - xyz_j2000_ts.time[i-1]

            xyz_dot = (xyz - xyz_prev) / dt

        Z_LOF = - xyz/np.linalg.norm(xyz)
        Y_LOF = - np.cross(xyz,xyz_dot)/np.linalg.norm(np.cross(xyz,xyz_dot))
        X_LOF = np.cross(Y_LOF,Z_LOF)

        R_J2000_LOF = np.array([X_LOF, Y_LOF, Z_LOF]).T

        # check
        u = np.dot(R_J2000_LOF,[0,0,1])

        #
        q_j2000_lof = Math.rotation_matrix_to_quat(R_J2000_LOF)

        q_j2000_lof_ts.append(xyz_j2000_ts.time[i], q_j2000_lof)

        q_j2000_sc = Math.quat_mult(q_j2000_lof, q_lof_sc)

        q_j2000_sc_ts.append(xyz_j2000_ts.time[i], q_j2000_sc)

    return q_j2000_lof_ts, q_j2000_sc_ts

def propagatePositionFromKeplerianElements(kepler_parameters,DT,dt, time_offset = 0, dateRef='2025-03-21T12:00:00Z'):

    '''
    This function propagates a spacecraft orbit defined with keplerian parameters and returns cartesian coordinates

    The initial position is represented by the kepler paratemers at a reference time
    The time offset parameter updates the mean anomaly from the reference time to the new starting time (t=0 in the timeseries)
    '''

    ## Keplerian parameters
    sma = kepler_parameters[0]
    ecc = kepler_parameters[1]
    inc = kepler_parameters[2]
    raan = kepler_parameters[3]
    aop = kepler_parameters[4]
    ta = kepler_parameters[5]

    ## Initialize timeseries

    xyz_j2000_ts = Timeseries('pos_j2000',dateRef)
    xyz_ecef_ts = Timeseries('pos_ecef',dateRef)
    ta_ts  = Timeseries('true_anomaly', dateRef)
    
    Nsteps = int(np.ceil(DT/dt));

    ## Orbital period
    T = orbitalPeriod(sma);

    ## Initial mean anomaly
    E0 = frames.eccentricAnomalyFromTrueAnomaly(ecc,ta);
    M0 = frames.meanAnomalyFromEccentricAnomaly(E0,ecc);

    for i in range(Nsteps):
        t = i*dt

        # Update mean anomaly
        M = np.mod(M0 + 2*np.pi*(t+time_offset)/T, 2*np.pi)
        ea = frames.eccentricAnomalyFromMeanAnomaly(M,ecc)
        ta = frames.trueAnomalyFromEccentricAnomaly(ea,ecc)

        ta_ts.append(t,ta);

        xyz = frames.kepler2cartesian(sma, ecc, inc, raan, aop, ta);
        xyz_j2000_ts.append(t,xyz)

        xyz_ecef_ts.append(t,frames.j2000_to_ecef(dateRef, xyz, time_offset+t))

    return xyz_j2000_ts, xyz_ecef_ts, ta_ts


class SpacecraftPosition(object):

    refTime = 'YYYY-MM-DDTHH:MM:SSZ'
    kepler_parameters = np.zeros(6)

    def __init__(self, refTime = '2025-03-21T12:00:00Z', kepler_parameters = np.zeros(6)):
        self.refTime = refTime
        self.kepler_parameters = kepler_parameters

    def propagate_from_start_date(self,startTime,DT,dt):
        # Time offset
        time_offset = Dates.time_offset(startTime, self.refTime)

        # propagate
        xyz_j2000_ts, xyz_ts,ta_ts = propagatePositionFromKeplerianElements(self.kepler_parameters,DT,dt, time_offset, self.refTime)

        xyz_j2000_ts.refTime = startTime
        xyz_ts.refTime = startTime
        ta_ts.refTime = startTime

        return xyz_j2000_ts,xyz_ts, ta_ts



class Spacecraft(object):
    """Spacecraft class"""

    id = ""
    norad_id = ""

    referencePosition = SpacecraftPosition()
    
    ta = Timeseries('true_anomaly');
    xyz = Timeseries('pos_ecef');
    xyz_j2000 = Timeseries('pos_j2000');
    h_long_lat = Timeseries('h_long_lat');

    q_j2000_lof = Timeseries('qj2000_lof');
    q_j2000_sc = Timeseries('q_j2000_sc');

    antenna_aperture = 0;
    
    def __init__(self, id, kepler_parameters=np.zeros(6), refTime = '2025-03-21T12:00:00Z'):
        """Constructor"""
        self.id = id
        self.referencePosition = SpacecraftPosition(refTime = refTime,kepler_parameters=kepler_parameters)
        
    def get_orbitalPeriod(self):
        return orbitalPeriod(self.referencePosition.kepler_parameters[0])

    def propagate(self,DT,dt, refDate = ""):

        if refDate == "":
            refDate = self.referencePosition.refTime

        self.xyz_j2000, self.xyz, self.ta = propagatePositionFromKeplerianElements(self.referencePosition.kepler_parameters,DT,dt, 0, refDate)

    def propagate_from_start_date(self, startDate, DT = 3600, dt = 60):
        self.xyz_j2000, self.xyz, self.ta = self.referencePosition.propagate_from_start_date(startDate,DT,dt)

    def set_default_attitude_pointing(self, q_lof_sc = [0,0,0,1]):

        self.q_j2000_lof, self.q_j2000_sc = compute_lof_sc_guidance_with_constant_bias(self.xyz_j2000, q_lof_sc)


    def get_position_h_long_lat(self):

        h_long_lat_ts = Timeseries('h_long_lat')
        h_long_lat_ts.refTime = self.xyz.refTime
        
        refTime = Time(self.xyz.refTime, scale='utc')

        # Calculate 
        for t,xyz in zip(self.xyz.time, self.xyz.data):

            time = refTime + t*u.s

            r,long,lat = frames.cartesian2spherical(xyz)

            h_long_lat_ts.append(t,[r-EARTH_RADIUS, long, lat])

        self.h_long_lat = h_long_lat_ts
        return h_long_lat_ts
