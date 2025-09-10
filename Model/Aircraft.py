import numpy as np

from MissionAnalysis.CommonTools.Timeseries import Timeseries

from datetime import datetime

## AERONAUTICAL UNITS

FEET_TO_KM = 0.30480*1e-3 # km

FL_TO_KM = 100*FEET_TO_KM # km
KM_TO_FL = 1/FL_TO_KM

KNOTS_TO_KPH = 1.852 # km/h

##


EARTH_R = 6378 # km

#

SECS_TO_H = 1/3600
H_TO_SECS = 3600

M_TO_KM = 1E-3


def propagate_position(pos_h_long_lat, vel_vx_vyz_track, DT=3600, dt=60):

    trajectory_h_long_lat = Timeseries()

    dt = dt*np.sign(DT)

    Nsteps = int(np.ceil(DT/dt))

    # initialize
    h = pos_h_long_lat[0]
    long = pos_h_long_lat[1]
    lat = pos_h_long_lat[2]

    # update speed - change units
    vx = vel_vx_vyz_track[0]/H_TO_SECS
    vyz = vel_vx_vyz_track[1]/H_TO_SECS
    track = vel_vx_vyz_track[2]


    # loop
    for i in range(Nsteps):

        t = i*dt
             
        # update h, long, lat derivatives
        r = EARTH_R + h
        h_dot = vx
        long_dot = vyz*np.sin(track)/r/np.cos(lat)
        lat_dot = vyz*np.cos(track)/r

        # update
        h = h + dt*h_dot
        long = long + dt*long_dot
        lat = lat + dt*lat_dot

        # store timeseries
        trajectory_h_long_lat.append(t, [h,long,lat])

    return trajectory_h_long_lat


class AircraftPosition(object):
    '''
    AircraftPosition object
    '''

    alt__km = 0
    long = 0
    lat = 0

    timestamp = 'YYYY-MM-DDTHH:MM:SSZ'  # UTC ISO 8601 format
    source = ''
    
    gspeed = 0  #km/h
    vspeed = 0  #km/h 
    track = 0   # True track (over ground) expressed in integer degrees as 0-360 taken clockwise from North

    def __init__(self, alt__km=0, long=0, lat=0, timestamp = '2000-01-01T11:59:07Z', source = 'unknown'):
        self.alt__km = alt__km
        self.long = long
        self.lat  = lat
        self.timestamp = timestamp
        self.source = source
        

    def set_trajectory(self, gspeed, vspeed, track):
        self.gspeed = gspeed
        self.vspeed = vspeed
        self.track  = track


    def propagate_position_from_start_date(self, startTime, DT = 3600, dt = 60):

        newTime = datetime.strptime(startTime, '%Y-%m-%dT%H:%M:%SZ')

        # Time offset
        if self.timestamp!='':
            timeRef = datetime.strptime(self.timestamp, '%Y-%m-%dT%H:%M:%SZ')
            time_offset = (newTime-timeRef).total_seconds()
        else: 
            self.timestamp = startTime
            time_offset = 0

        # Refresh start position
        h_long_lat_start = [self.alt__km, self.long, self.lat]
        if abs(time_offset)>=dt:
            h_long_lat_ts = propagate_position(h_long_lat_start,
                                           [self.vspeed, self.gspeed, self.track], 
                                           DT=time_offset, dt=60)
            toffset, h_long_lat_start = h_long_lat_ts.get_data_at_index(-1)
        
        # Propagate from new start position
        h_long_lat_ts = propagate_position(h_long_lat_start,
                                       [self.vspeed, self.gspeed, self.track], 
                                       DT=DT, dt=dt)
        h_long_lat_ts.refTime = startTime

        return h_long_lat_ts






class Aircraft(object):
    '''
    Aircraft Model
    '''

    icao_id = ''
    flight_id = ''
    airline = ''

    departure_airport  = ''
    arrival_airport    = ''

    reference_position = AircraftPosition()

    trajectory = Timeseries()

    def __init__(self, icao_id):
        self.icao_id = icao_id

    def set_reference_position(self, alt__km, long, lat, timestamp = '2000-01-01T11:59:07Z', source = 'unknown'):
        self.reference_position = AircraftPosition(alt__km, long, lat, timestamp, source)

    def set_reference_trajectory(self, gspeed, vspeed, track):
        self.reference_position.set_trajectory(gspeed, vspeed, track)

    def get_position_at_time(self, time):
        return self.reference_position.get_position_at_time(time)

    def propagate_position(self, startTime, DT = 3600, dt = 60):
        self.trajectory = self.reference_position.propagate_position_from_start_date(startTime, DT, dt)
        return self.trajectory

    def return_as_dict(self, time, h_long_lat):

        gbs_flag = '---'

        if h_long_lat[0] == 0:
            gbs_flag = 'GBS'

        aircraft = {
            'icao_id': self.icao_id, 
            'flight_id': self.flight_id, 
            'time': time,
            'gbs': gbs_flag,
            'altitude (km)': h_long_lat[0],
            'long (deg)': h_long_lat[1],
            'lat (deg)': h_long_lat[2]}

        return aircraft


    def __str__(self):

        return f'Aircraft(id = {self.icao_id}, position = {self.reference_position}'





