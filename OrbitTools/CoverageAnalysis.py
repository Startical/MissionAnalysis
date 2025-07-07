import numpy as np
import pandas as pd

from CommonTools.Constants import EARTH_RADIUS
from OrbitTools.FrameTransformations import FrameTransformations as frames
from datetime import datetime, timezone, timedelta


def antenna_swath(h,H, beta):
    '''
    This function calculates the antenna swath (footprint) at any altitude (H, H=0 = sea level)
    given the terminal altitude (h), the antenna aperture halfcone angle (beta)

    INPUTS: 
    - h     : terminal altitude (km)
    - H     : target altitude (km)
    - beta  : antenna aperture half cone angle (rad)

    OUTPUTS:
    - alpha : footprint angle at H (rad)
    - c     : footprint radius at H (km)
    '''

    # Calculate limit of Earth horizon, betalim
    hh = (EARTH_RADIUS+h)/(EARTH_RADIUS+H);
    betalim = np.arcsin(1/hh)
    
    if beta <= betalim:

        alpha = 0;
        while np.abs(np.sin(alpha)-(hh-np.cos(alpha))*np.tan(beta)) > 1e-6:
            alpha = np.arcsin((hh-np.cos(alpha))*np.tan(beta));

    else: 
        betamin = np.arcsin(EARTH_RADIUS/(EARTH_RADIUS+h))

        # First solution, before tangent point
        #alpha1 = 0;
        #
        #while np.abs(np.sin(alpha1)-(hh-np.cos(alpha1))*np.tan(betamin)) > 1e-6:
        #    alpha1 = np.arcsin((hh-np.cos(alpha1))*np.tan(betamin));

        # Second solution, after tangent point
        deltaAlpha = np.pi/2/1000
        alpha2 = np.pi/2-betalim;
        while (np.sin(alpha2)-(hh-np.cos(alpha2))*np.tan(betamin)) > 0:
            alpha2 = alpha2+deltaAlpha

        alpha = alpha2
                
    c = (EARTH_RADIUS+H)*alpha

    return alpha, c



def check_visibility_lower_upper_terminal(lower_terminal_pos, upper_terminal_pos, lower_terminal_fov=np.pi, upper_terminal_fov=np.pi):

    """
    Check if the lower terminal is visible from the upper terminal and vice versa.
    The function returns True if the lower terminal is visible from the upper terminal
    and False otherwise.
    terminal_pos = [h, long, lat]
    """
    
    is_visible = False

    r_up = upper_terminal_pos[0]+EARTH_RADIUS
    ecef_xyz_up = frames.spherical2cartesian(r_up,upper_terminal_pos[1],upper_terminal_pos[2]) 
    ecef_u_up = ecef_xyz_up/r_up

    r_lw = lower_terminal_pos[0]+EARTH_RADIUS
    ecef_xyz_lw = frames.spherical2cartesian(r_lw,lower_terminal_pos[1],lower_terminal_pos[2])
    ecef_u_lw = ecef_xyz_lw/r_lw

    ecef_xyz_lw_up = ecef_xyz_up - ecef_xyz_lw
    ecef_u_lw_up = ecef_xyz_lw_up/np.linalg.norm(ecef_xyz_lw_up)

    # horizon limit
    gamma_lim = np.arcsin(np.clip(EARTH_RADIUS/r_lw,-1.0,+1.0))
    # check if upper terminal is above the Earth's horizon
    if (angle_u_v(-ecef_u_lw, ecef_u_lw_up) > gamma_lim):
        # satellite is above the horizon, check if visible within antenna apertures
        if (angle_u_v(ecef_u_lw, ecef_u_lw_up)< lower_terminal_fov) and (angle_u_v(-ecef_u_up, -ecef_u_lw_up)< upper_terminal_fov):
            is_visible = True
            
    return is_visible


def check_ground_station_spacecraft_visibility(spacecraft, ground_stations):

    refTime = spacecraft.h_long_lat.refTime

    visibility_timelapse = [[] for _ in range(len(ground_stations))]

    for t, sc_h_long_lat in zip(spacecraft.h_long_lat.time, spacecraft.h_long_lat.data):

        for idx, ground_station in enumerate(ground_stations):
            if check_visibility_lower_upper_terminal(ground_station.get_position(), sc_h_long_lat, lower_terminal_fov=ground_station.antenna_aperture, upper_terminal_fov=np.pi):
                visibility_timelapse[idx].append(t)


    # Create result structure
    
    gs_visibility = pd.DataFrame(columns=['Ground Station', 'Visibility opportunity', 'Visibility starts', 'Visibility ends', 'Duration (min)'])
    rows = []

    for idx, ground_station in enumerate(ground_stations):

        is_visible = False
        t_start = ''
        t_end = ''
        duration = 0

        if len(visibility_timelapse[idx])>0:
            is_visible = True
            t_start =  (datetime.strptime(refTime, '%Y-%m-%dT%H:%M:%SZ') + timedelta(seconds=visibility_timelapse[idx][0])).strftime('%Y-%m-%dT%H:%M:%SZ')
            t_end =  (datetime.strptime(refTime, '%Y-%m-%dT%H:%M:%SZ') + timedelta(seconds=visibility_timelapse[idx][-1])).strftime('%Y-%m-%dT%H:%M:%SZ')
            duration = visibility_timelapse[idx][-1] - visibility_timelapse[idx][0]

        new_row = {
            'Ground Station': ground_station.gs_id,
            'Visibility opportunity': is_visible,
            'Visibility starts': t_start,
            'Visibility ends': t_end,
            'Duration (min)': duration/60
        }

        rows.append(new_row)

    gs_visibility = pd.concat([gs_visibility, pd.DataFrame(rows)], ignore_index=True)

    return gs_visibility, visibility_timelapse

def angle_u_v(u,v):
    u_unit = u/np.linalg.norm(u)
    v_unit = v/np.linalg.norm(v)

    f = np.dot(u_unit,v_unit)

    if f<-1.0:
        f=-1.0
    if f>1.0:
        f=1.0

    return np.arccos(f)




