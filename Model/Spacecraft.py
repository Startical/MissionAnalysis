import numpy as np
from OrbitTools.Timeseries import Timeseries
from OrbitTools.Constants import EARTH_MU
from OrbitTools.FrameTransformations import FrameTransformations as frames

class Spacecraft(object):
    """Spacecraft class"""

    id = ""
    kepler_parameters = np.zeros(6);
    T = 0; # orbital period
    ta = Timeseries();
    xyz = Timeseries();

    antenna_aperture = 0;
    
    def __init__(self, id, kepler_parameters):
        """Constructor"""
        self.id = id
        self.kepler_parameters = kepler_parameters

    def get_orbitalPeriod(self):
        """Get the orbital period"""

        sma = self.kepler_parameters[0]

        self.T = 2*np.pi*np.sqrt(sma**3/EARTH_MU)

        return self.T


    def propagate(self,t,dt):

        self.xyz = Timeseries()

        Nsteps = int(np.ceil(t/dt));

        ## Orbital period
        T = self.get_orbitalPeriod();

        ## Initial mean anomaly
        E0 = frames.eccentricAnomalyFromTrueAnomaly(self.kepler_parameters[1],self.kepler_parameters[5]);
        M0 = frames.meanAnomalyFromEccentricAnomaly(E0,self.kepler_parameters[1]);

        for i in range(Nsteps):
            t = i*dt

            # Update mean anomaly
            M = M0 + 2*np.pi*t/T
            ea = frames.eccentricAnomalyFromMeanAnomaly(M,self.kepler_parameters[1])
            ta = frames.trueAnomalyFromEccentricAnomaly(ea,self.kepler_parameters[1])

            self.ta.append(t,ta);

            xyz = frames.kepler2cartesian(self.kepler_parameters[0],self.kepler_parameters[1],self.kepler_parameters[2],self.kepler_parameters[3],self.kepler_parameters[4],ta);
            self.xyz.append(t,xyz)

        return self.xyz