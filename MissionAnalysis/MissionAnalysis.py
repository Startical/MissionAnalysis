

import numpy as np
from ConstellationDesign import Constellation
import matplotlib
matplotlib.use('Qt5Agg')  # or 'Qt5Agg' if you prefer

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Circle

import cartopy.crs as ccrs
import cartopy.feature as cfeature

from OrbitTools.FrameTransformations import FrameTransformations as frames
from OrbitTools.Constants import EARTH_RADIUS

def plot_earth(constellation):
    """Plot the Earth in longitude/latitude coordinates"""
    fig, ax = plt.subplots(figsize=(10, 5), subplot_kw={'projection': ccrs.PlateCarree()})

    # Add coastlines and borders
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle=':')

    # Add gridlines
    gl = ax.gridlines(draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False

    # Define a colormap
    colormap = cm.get_cmap('viridis', constellation.m)

    # Plot satellites
    for sc in constellation.spacecraft:
        xyz = sc.xyz.get_data_at_time(0)
        R,long,lat = frames.cartesian2spherical(xyz)
        i_index = int(sc.id.split('_(')[1].split(',')[0])  # Extract the i index from the spacecraft ID
        color = colormap(i_index / constellation.m)  # Get the color from the colormap
        ax.scatter(long*180/np.pi, lat*180/np.pi, color=color, s=10, transform=ccrs.PlateCarree())
        # Add a disk (circle) centered on each point
        r = constellation.antenna_swath()*180/np.pi;
        circle = Circle((long*180/np.pi, lat*180/np.pi), radius=r, color=color, alpha=0.2, transform=ccrs.PlateCarree())
        ax.add_patch(circle)

    plt.title("Earth in Longitude/Latitude Coordinates")
    plt.show()

if __name__ == "__main__":

    print("hello");

    constellation = Constellation("Default")

    constellation.h = 650
    constellation.inc = 80*np.pi/180
    constellation.n = 16
    constellation.m = 9
    constellation.walker_option = "star"

    constellation.H = 18; # FL600

    constellation.antenna_aperture = np.asin((EARTH_RADIUS+constellation.h)/(EARTH_RADIUS+constellation.H));

    constellation.initializeSpacecraft()

    plot_earth(constellation)

    print("hello")



