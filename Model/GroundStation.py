
import numpy as np
from OrbitTools.Constants import EARTH_RADIUS


class GroundStation(object):
    '''
    Ground station object
    '''

    gs_id = ''
    h = 0
    long = 0
    lat = 0

    services = ''

    antenna_aperture = np.deg2rad(85) # default value, 5º above horizon

    def __init__(self, gs_id='', alt = 0, long=0, lat=0, services='', antenna_aperture = np.deg2rad(85)):
        self.gs_id = gs_id
        self.h = alt
        self.long = long
        self.lat = lat
        self.services = services
        self.antenna_aperture = antenna_aperture

    def get_position(self):
        return [self.h, self.long, self.lat]
    
    def get_great_circle_angle_at_h(self, h_target):
        ''' 
        Calculates great circle angle from ground station to target at altitude h considering maximum antenna aperture
        '''
        el = np.pi/2 - self.antenna_aperture
        return self.get_great_circle_angle_to_target(h_target, el)

    def get_great_circle_angle_to_target(self, h_target, el):
        ''' 
        Calculates great circle angle from ground station to target at altitude h and elevation el
        '''

        beta_target = np.arcsin((EARTH_RADIUS + self.h)/(EARTH_RADIUS + h_target)*np.cos(el))
        theta = np.pi/2 - beta_target - el

        return theta




