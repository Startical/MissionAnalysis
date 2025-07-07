import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Polygon
from geopy.distance import distance
from geopy.point import Point
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from MissionAnalysis.OrbitTools.Constants import *
import numpy as np

def initialize_Earth_2D_plot():

    ## Print Earth Map
    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
    
    # Add coastlines and borders
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    
    # Add gridlines
    gl = ax.gridlines(draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False

    return fig, ax


def plot_angle_variable(time,angles):
    '''
    This function inserts a nan value in an angle series to create a discontinuity for plots
    '''
    angles_with_nan = []
    time_with_nan = []
    for i in range(len(angles) - 1):
        angles_with_nan.append(angles[i])
        time_with_nan.append(time[i])
        if abs(angles[i] - angles[i + 1]) > 180:
            angles_with_nan.append(np.nan)
            time_with_nan.append(np.nan)
    angles_with_nan.append(angles[-1])
    time_with_nan.append(time[-1])

    return time_with_nan, angles_with_nan

def haversine(lat1, lon1, lat2, lon2):

    lat1, lon1 = np.radians(lat1), np.radians(lon1)
    lat2, lon2 = np.radians(lat2), np.radians(lon2)

    dlat = lat2 - lat1
    dlon = lon2 - lon1

    a = np.sin(dlat/2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.0)**2
    c = 2 * np.arcsin(np.sqrt(a))
    return EARTH_RADIUS * c

def circle_projection(long, lat, az,radius):
 
    lon_new = []
    lat_new = []
    for theta in az:
        lon_aux = []
        lat_aux = []
        for r in radius:
            coords = distance(kilometers=r*EARTH_RADIUS).destination(point=Point(lat,long), bearing=np.rad2deg(theta))
            if np.deg2rad(coords.latitude) > 90:
                print("Warning")
            lon_aux.append(np.deg2rad(coords.longitude))
            lat_aux.append(np.deg2rad(coords.latitude))
        
        lon_new.append(lon_aux)
        lat_new.append(lat_aux)
            
    return np.array(lon_new),np.array(lat_new)

def plot_circle_projection(ax, long, lat, r, color='b', alpha = 0.2):

    r1 = r
    r2 = r

    N = 18
    vertices = []

    for i in range(N):
        theta = 2*np.pi*i/N
        lat1 = lat + r*np.sin(theta)
        dlong = r*np.cos(theta)

        if lat1 > 180-(lat+r):
            lat1 = 180-(lat+r)
            dlong = (dlong+180+180) % 360 - 180
        elif lat1 < -180-(lat-r):
            lat1 = -180-(lat-r)
            dlong = (dlong+180+180) % 360 - 180


        dlong = dlong/np.cos(np.deg2rad(lat1))

        long1 = np.clip(long+dlong, -180, 180) #(long + dlong +180) % 360 -180

        # 
        vertices.append([long1,lat1])

        r1 = min(max(r1,dlong),180)

    if (lat+r)>90 or (lat-r)<-90:
        # sort values and add corners to ensure that the area containing the centre (long,lat) is printed by the Polygon
        vertices = sorted(vertices, key=lambda x: x[0])
        vertices.insert(0,[-180, np.sign(lat)*90])
        vertices.append([180, np.sign(lat)*90])

    projected_circle = Polygon(vertices, color=color, alpha=alpha, transform=ccrs.PlateCarree())

    ax.add_patch(projected_circle)

    return r1,r2

def plot_ground_station(ax, gs, theta, gs_color = [0,1,0]):
    
    # Print the ground station coverage
    r1,r2 = plot_circle_projection(ax, np.rad2deg(gs.long), np.rad2deg(gs.lat), np.rad2deg(theta), color=gs_color, alpha=0.2)
    # Print the ground station position
    ax.scatter(np.rad2deg(gs.long), np.rad2deg(gs.lat), marker='d', color = gs_color, label=gs.gs_id, transform=ccrs.PlateCarree())

    return r1, r2

def initialize_ground_track_plot(spacecraft):
    
    fig, ax = initialize_Earth_2D_plot()

    ## Print satellite trajectory
    sc_h_long_lat_ts = spacecraft.get_position_h_long_lat()

    lat, long = plot_angle_variable(np.rad2deg(sc_h_long_lat_ts.data[:,2]), np.rad2deg(sc_h_long_lat_ts.data[:,1]))

    ax.plot(long, lat, 'k--', label=spacecraft.id)

    for idx in np.arange(0,len(sc_h_long_lat_ts.data),10):
        arrow_length = 5
        dlong = sc_h_long_lat_ts.data[idx+1,1]-sc_h_long_lat_ts.data[idx,1]
        dlat = sc_h_long_lat_ts.data[idx+1,2]-sc_h_long_lat_ts.data[idx,2]
        ax.arrow(np.rad2deg(sc_h_long_lat_ts.data[idx,1]),np.rad2deg(sc_h_long_lat_ts.data[idx,2]),
                dlong/np.linalg.norm([dlong,dlat]) * arrow_length, dlat/np.linalg.norm([dlong,dlat])  * arrow_length, color='k', head_width=1, head_length=2, transform=ccrs.PlateCarree())
    
    ax.set_title(f'ground_track_{spacecraft.id}\n{spacecraft.xyz.refTime}')
    fig.name = f'ground_track_{spacecraft.id}'

    return fig, ax