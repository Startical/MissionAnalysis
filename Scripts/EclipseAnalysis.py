
import numpy as np

from Model.SpaceEnvironment import sun_direction_inertial
from CommonTools.Constants import EARTH_RADIUS, EARTH_MU
from OrbitTools.FrameTransformations import FrameTransformations as frames

import matplotlib.pyplot as plt
import matplotlib.cm as cm


def calculate_eclipse_duration(RAAN, inc, month, h=650):
    '''
    Calculates the eclipse duration in minutes
    '''
    
    # sun elevation over orbital plane
    beta = np.abs(sun_elevation_orbital_plane(RAAN, inc, month))

    # calculate max elevation for eclipse
    beta_lim = np.arcsin(EARTH_RADIUS / (EARTH_RADIUS + h))

    if beta > beta_lim:
        return 0 # no eclipse

    # orbital period
    T = 2 * np.pi * np.sqrt((EARTH_RADIUS + h)**3 / EARTH_MU)

    # calculate arc segment with eclipse
    arc_segment = 2 * np.arccos(np.clip(np.cos(beta_lim) / np.cos(beta),-1.0,+1.0))
    
    # calculate eclipse duration in minutes
    eclipse_duration = (arc_segment / (2 * np.pi)) * T / 60

    return eclipse_duration

def sun_elevation_orbital_plane(RAAN, inc, month):
    '''
    Calculates sun elevation angle wrt the orbital plane (beta angle)
    '''

    # Calculate the sun direction in ECI frame
    sun_direction = sun_direction_inertial(month)

    # Calculate the orbital plane normal vector
    orbit_N = np.dot(frames.orbitPlane2Cartesian(inc,RAAN,0),[0,0,1])

    # Calculate the sun elevation angle
    beta = np.arcsin(np.clip(np.dot(sun_direction, orbit_N),-1.0,+1.0))

    return beta

def plot_sun_elevation_orbital_plane():

    # Define parameters

    inc = np.array([45,60,75,90])*np.pi/180

    RAAN = np.linspace(0,2*np.pi,36)
    month = np.linspace(1,12,12)

    # Create a meshgrid for RAAN and month
    month_grid, RAAN_grid = np.meshgrid(month, RAAN)

    # Initialize figure
    fig,ax = plt.subplots(2,2)
    levels = np.linspace(0,90,19)

    # Loop on inclination
    for idx in range(4):

        beta = np.zeros(month_grid.shape)
        for i in range(len(month)):    
            for j in range(len(RAAN)):
            
                beta[j,i] = sun_elevation_orbital_plane(RAAN[j], inc[idx], month[i])
                
        # Plot
        contour = ax[idx//2, idx%2].contourf(month_grid, RAAN_grid*180/np.pi, np.abs(beta)*180/np.pi, levels=levels, cmap='jet')
        ax[idx//2, idx%2].set_title(f'Inc {inc[idx]*180/np.pi:.0f}deg')
        ax[idx//2, idx%2].set_ylabel('RAAN (deg)')
        ax[idx//2, idx%2].set_xlabel('Month')
        plt.colorbar(contour, ax=ax[idx//2, idx%2], orientation='horizontal')

    return fig


def plot_eclipse_duration(h, inc):

    # Define parameters
    RAAN = np.linspace(0,2*np.pi,36)
    month = np.linspace(1,12,12)

    # Create a meshgrid for RAAN and month
    month_grid, RAAN_grid = np.meshgrid(month, RAAN)

    # Initialize figure
    fig,ax = plt.subplots()

    eclipse_duration = np.zeros(month_grid.shape)
    for i in range(len(month)):    
        for j in range(len(RAAN)):
            
            eclipse_duration[j,i] = calculate_eclipse_duration(RAAN[j], inc, month[i], h)
                
    # Plot
    contour = ax.contourf(month_grid, RAAN_grid*180/np.pi, eclipse_duration, cmap='jet')
    ax.set_title(f'h = {h}km; Inc = {inc*180/np.pi:.0f}deg')
    ax.set_ylabel('RAAN (deg)')
    ax.set_xlabel('Month')
    cbar = plt.colorbar(contour, ax=ax, orientation='vertical')
    cbar.set_label('Eclipse duration (min)')


    return fig


def plot_eclipse_duration_year(h, inc, RAAN):

    # Define parameters
    month = np.linspace(1,12,365)


    # Initialize figure
    fig,ax = plt.subplots()

    eclipse_duration = np.zeros(month.shape)
    for i in range(len(month)):                
        eclipse_duration[i] = calculate_eclipse_duration(RAAN, inc, month[i], h)
                
    # Plot
    ax.plot(month, eclipse_duration)
    ax.set_title(f'h = {h}km; Inc = {inc*180/np.pi:.0f}deg; RAAN = {RAAN*180/np.pi:.0f}deg')
    ax.set_ylabel('Eclipse duration (min)')
    ax.set_xlabel('Month')
    ax.grid()

    return fig

if __name__ == "__main__":
    fig1 = plot_sun_elevation_orbital_plane()
    fig2 = plot_eclipse_duration(650, 75*np.pi/180)

    # Eclipse duration in GEO
    fig0 = plot_eclipse_duration_year(36000, 0, 0)

    plt.show()



