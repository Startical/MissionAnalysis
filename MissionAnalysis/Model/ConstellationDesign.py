from encodings import mac_farsi
from Model.Spacecraft import Spacecraft
from OrbitTools.Constants import EARTH_RADIUS
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')  # or 'Qt5Agg' if you prefer

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Circle
from matplotlib.colors import ListedColormap, BoundaryNorm, LinearSegmentedColormap

import cartopy.crs as ccrs
import cartopy.feature as cfeature

from OrbitTools.FrameTransformations import FrameTransformations as frames


class Constellation(object):
    """Constellation class"""

    id = ""
    spacecraft = []

    H = 0; # flight level
    h = 650;
    inc = 90*np.pi/180;
    walker_option = "star";
    n = 16; # number of satellites per plane
    m = 9; # number of orbital planes
    N = n*m;

    antenna_aperture = 60*np.pi/180;

    def __init__(self, id):
        self.id = id;

    def __init__(self, id, h, inc, m, n, walker_option, antenna_aperture, H):
        self.id = id
        self.h = h
        self.inc = inc
        self.m = m
        self.n = n
        self.walker_option = walker_option
        self.antenna_aperture = antenna_aperture
        self.H = H
        self.N = m*n

    def idFullContext(self):
        return f'{self.id}_{self.h:.0f}_{self.inc*180/np.pi:.0f}_{self.m}_{self.n}_{self.walker_option}_{self.antenna_aperture*180/np.pi:.0f}_FL{self.H/3*100:.0f}'

    def antenna_swath(self):

        hh = (EARTH_RADIUS+self.h)/(EARTH_RADIUS+self.H);
        gamma = np.acos(1/hh);

        while np.sin(gamma)-(hh-np.cos(gamma))*np.tan(self.antenna_aperture) > 1e-6:
            gamma = np.asin((hh-np.cos(gamma))*np.tan(self.antenna_aperture));

        return gamma


    def planesRAANdistribution(self):
        raanDistribution = 2*np.pi
        if self.walker_option == "star":
            raanDistribution = np.pi

        return raanDistribution

    def initializeSpacecraft(self):
        
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
                newSC.propagate(newSC.get_orbitalPeriod(),newSC.get_orbitalPeriod()/20)
                self.spacecraft.append(newSC)

    def plot_constellation_groundTrack(self):
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
        colormap = cm.get_cmap('viridis', self.m)
    
        # Plot satellites
        for sc in self.spacecraft:
            xyz = sc.xyz.get_data_at_time(0)
            dt, xyz1 = sc.xyz.get_data_at_index(1)
    
            R,long,lat = frames.cartesian2spherical(xyz)
            R1, long1, lat1 = frames.cartesian2spherical(xyz1)
    
            i_index = int(sc.id.split('_(')[1].split(',')[0])  # Extract the i index from the spacecraft ID
            color = colormap(i_index / self.m)  # Get the color from the colormap
            
            ax.scatter(long*180/np.pi, lat*180/np.pi, color=color, s=10, transform=ccrs.PlateCarree())
            ax.text(long*180/np.pi, lat*180/np.pi, sc.id, fontsize=4, verticalalignment='bottom')
            arrow_length = 5
            dlong = long1-long
            dlat = lat1-lat
            ax.arrow(long * 180 / np.pi, lat * 180 / np.pi, dlong/np.linalg.norm([dlong,dlat]) * arrow_length, dlat/np.linalg.norm([dlong,dlat])  * arrow_length, color=color, head_width=1, head_length=2, transform=ccrs.PlateCarree())
    
            # Add a disk (circle) centered on each point
            r = self.antenna_swath()*180/np.pi;
            circle = Circle((long*180/np.pi, lat*180/np.pi), radius=r, color=color, alpha=0.2, transform=ccrs.PlateCarree())
            ax.add_patch(circle)
    
        plt.title(f'Constellation coverage: {self.h} km, {self.inc*180/np.pi:.0f} deg, {self.m} x {self.n} , walker {self.walker_option}')
        
        fig.name = f'{self.idFullContext()}_fig1'

        return fig;
    
    def plot_constellation_coverage(self):
        H = self.H
    
        delta = np.pi/50
        long = np.arange(0,2*np.pi+delta, delta)
        lat = np.arange(-np.pi/2, np.pi/2+delta, delta)
    
        long_grid, lat_grid = np.meshgrid(long, lat)
        n = np.zeros_like(long_grid)
    
        for i in range(len(long)): 
            for j in range(len(lat)):
    
                xyz_P = frames.spherical2cartesian(EARTH_RADIUS + H, long[i], lat[j]);
                R = np.linalg.norm(xyz_P)
                u_P = xyz_P/R
                gamma_lim = np.asin(EARTH_RADIUS/R)
    
                for sc in self.spacecraft:
    
                    xyz_SAT = sc.xyz.get_data_at_time(0)
                    u_SAT = xyz_SAT/np.linalg.norm(xyz_SAT)
                    xyz_P_SAT = xyz_SAT - xyz_P
                    u_P_SAT = xyz_P_SAT/np.linalg.norm(xyz_P_SAT)
    
                    # check if satellite is above the horizon
                    if (angle_u_v(-u_P, u_P_SAT) > gamma_lim):
                        # satellite is above the horizon, check if visible within antenna aperture
                        if (angle_u_v(-u_SAT, -u_P_SAT)< self.antenna_aperture):
                            n[j,i] = n[j,i]+1
        
        fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
    
        # Define contour levels
        levels = np.arange(0,6,1)
        levels = np.append(levels, np.max(n))
    
        # Create a custom colormap
        coolwarm = plt.get_cmap('coolwarm')
        colors = [(1, 1, 1)] + [coolwarm(i) for i in range(coolwarm.N)]
        cmap = ListedColormap(colors)
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    
        # Plot coverage as contour
        contour = ax.contourf(long_grid*180/np.pi, lat_grid*180/np.pi, n, levels=levels, cmap=cmap, norm=norm)
    
        # Add coastlines and borders
        ax.coastlines()
        ax.add_feature(cfeature.BORDERS, linestyle=':')
    
        # Add gridlines
        gl = ax.gridlines(draw_labels=True)
        gl.top_labels = False
        gl.right_labels = False
    
        # Add colorbar with title
        cbar = plt.colorbar(contour, ax=ax, orientation='horizontal', pad=0.05)
        cbar.set_label('Number of visible satellites')
    
        plt.title(f'Constellation coverage: {self.h} km, {self.inc*180/np.pi:.0f} deg, {self.m} x {self.n} , walker {self.walker_option}')
        fig.name = f'{self.idFullContext()}_fig2'

        return fig;

def angle_u_v(u,v):
    u_unit = u/np.linalg.norm(u)
    v_unit = v/np.linalg.norm(v)

    f = np.dot(u_unit,v_unit)

    if f<-1.0:
        f=-1.0
    if f>1.0:
        f=1.0

    return np.acos(f)

    
        

