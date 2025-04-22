
import numpy as np


def sun_direction_inertial(month):
    '''
    Simplified model to calculate the sun direction in Earth Centered Intertial frame
    '''

    # Declination
    eps = 23.439292 * np.pi / 180  # obliquity of the ecliptic in radians

    # The sun direction is calculated as a function of the month
    sun_angle = 2 * np.pi * (month - 3) / 12

    # Sun direction in ECI frame
    sun_direction = np.array([np.cos(sun_angle), np.sin(sun_angle) * np.cos(eps), np.sin(sun_angle) * np.sin(eps)])

    return sun_direction

if __name__ == "__main__":
    # Test the function
    month = 6  # June
    sun_direction = sun_direction_inertial(month)
    print("Sun direction in ECI frame:", sun_direction)