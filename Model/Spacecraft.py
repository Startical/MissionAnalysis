import numpy as np
from OrbitTools.Timeseries import Timeseries
from OrbitTools.Constants import EARTH_MU, EARTH_RADIUS, EARTH_W
from OrbitTools.FrameTransformations import FrameTransformations as frames
from datetime import datetime, timedelta
from astropy.coordinates import CartesianRepresentation, EarthLocation, ITRS, GCRS
from astropy.time import Time
import astropy.units as u



def orbitalPeriod(sma):
    """Compute the orbital period"""

    T = 2*np.pi*np.sqrt(sma**3/EARTH_MU)

    return T


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

    xyz_j2000_ts = Timeseries()
    xyz_ecef_ts = Timeseries()
    ta_ts  = Timeseries()
    
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
        timeRef = datetime.strptime(self.refTime, '%Y-%m-%dT%H:%M:%SZ')
        newTime = datetime.strptime(startTime, '%Y-%m-%dT%H:%M:%SZ')

        time_offset = (newTime-timeRef).total_seconds()

        # propagate
        xyz_j2000_ts, xyz_ts,ta_ts = propagatePositionFromKeplerianElements(self.kepler_parameters,DT,dt, time_offset, self.refTime)

        xyz_ts.refTime = startTime
        ta_ts.refTime = startTime

        return xyz_j2000_ts,xyz_ts, ta_ts



class Spacecraft(object):
    """Spacecraft class"""

    id = ""
    norad_id = ""

    referencePosition = SpacecraftPosition()
    
    ta = Timeseries();
    xyz = Timeseries();
    xyz_j2000 = Timeseries();
    h_long_lat = Timeseries();

    antenna_aperture = 0;
    
    def __init__(self, id, kepler_parameters=np.zeros(6), refTime = '2025-03-21T12:00:00Z'):
        """Constructor"""
        self.id = id
        self.referencePosition = SpacecraftPosition(refTime = refTime,kepler_parameters=kepler_parameters)
        
    def get_orbitalPeriod(self):
        return orbitalPeriod(self.referencePosition.kepler_parameters[0])

    def propagate(self,DT,dt, refDate = referencePosition.refTime):
        self.xyz_j2000, self.xyz, self.ta = propagatePositionFromKeplerianElements(self.referencePosition.kepler_parameters,DT,dt, 0, refDate)

    def propagate_from_start_date(self, startDate, DT = 3600, dt = 60):
        self.xyz_j2000, self.xyz, self.ta = self.referencePosition.propagate_from_start_date(startDate,DT,dt)

    def get_position_h_long_lat(self):

        h_long_lat_ts = Timeseries()
        h_long_lat_ts.refTime = self.xyz.refTime
        
        refTime = Time(self.xyz.refTime, scale='utc')

        # Calculate 
        for t,xyz in zip(self.xyz.time, self.xyz.data):

            time = refTime + t*u.s

            r,long,lat = frames.cartesian2spherical(xyz)

            h_long_lat_ts.append(t,[r-EARTH_RADIUS, long, lat])

        self.h_long_lat = h_long_lat_ts
        return h_long_lat_ts
