from Model.Spacecraft import Spacecraft
from CommonTools.Constants import EARTH_RADIUS

import numpy as np

class Constellation(object):
    """Constellation class"""

    name = ""
    spacecraft = []

    H = 0; # flight level
    h = 650;
    inc = 90*np.pi/180;
    walker_option = "star";
    n = 16; # number of satellites per plane
    m = 9; # number of orbital planes
    N = n*m;

    antenna_aperture = 60*np.pi/180;

    def __init__(self, name):
        self.name = name;

    def __init__(self, name, h, inc, m, n, walker_option, antenna_aperture, H):
        self.name = name
        self.h = h
        self.inc = inc
        self.m = m
        self.n = n
        self.walker_option = walker_option
        self.antenna_aperture = antenna_aperture
        self.H = H
        self.N = m*n

    def idFullContext(self):
        return f'{self.name}_{self.h:.0f}_{self.inc*180/np.pi:.0f}_{self.m}_{self.n}_{self.walker_option}_{self.antenna_aperture*180/np.pi:.0f}_FL{self.H/3*100:.0f}'

    def antenna_swath(self):

        hh = (EARTH_RADIUS+self.h)/(EARTH_RADIUS+self.H);
        gamma = np.arccos(1/hh);

        while np.sin(gamma)-(hh-np.cos(gamma))*np.tan(self.antenna_aperture) > 1e-6:
            gamma = np.arcsin((hh-np.cos(gamma))*np.tan(self.antenna_aperture));

        return gamma


    def planesRAANdistribution(self):
        raanDistribution = 2*np.pi
        if self.walker_option == "star":
            raanDistribution = np.pi

        return raanDistribution

    def initializeSpacecraft(self,DT,dt):
        
        # Initialize spacecraft parameters
        defaultKeplerParams = np.array([EARTH_RADIUS+self.h, 0, self.inc, 0, 0, 0])
        self.spacecraft = []

        for i in range(self.m):  # loop over planes
            for j in range(self.n):  # loop over satellites in each plane
                sc_id = f"SC_{i,j}"
                kepler_parameters = defaultKeplerParams

                kepler_parameters[3] = self.planesRAANdistribution()/self.m*i
                kepler_parameters[5] = 2*np.pi/self.n*(j+i/2) % (2 * np.pi)

                newSC = Spacecraft(sc_id, kepler_parameters)
                newSC.propagate(DT,dt)
                self.spacecraft.append(newSC)




    

