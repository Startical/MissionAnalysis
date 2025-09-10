import numpy as np

from datetime import datetime

from CommonTools.Constants import EARTH_W

from astropy.time import Time


class FrameTransformations(object):
    """This class performs frame transformations"""



    def xyz_coordinates_in_plane(sma,ecc,TA):
        """This method calculates the xyz coordinates in the orbital plane"""

        r = sma*(1 - ecc**2)/(1 + ecc*np.cos(TA))

        xyz = np.zeros(3); #initialize the array
        xyz[0] = r*np.cos(TA)
        xyz[1] = r*np.sin(TA)
        xyz[2] = 0
        return xyz

    def orbitPlane2Cartesian(inc, RAAN, AOP):
        """This method calculates the rotation matrix from the orbital plane frame (3) to cartesian (0)"""

        R_0_1 = np.array([
            [np.cos(RAAN), -np.sin(RAAN), 0],
            [np.sin(RAAN), np.cos(RAAN), 0],
            [0, 0, 1]
        ])
        R_1_2 = np.array([
            [1, 0, 0],
            [0, np.cos(inc), -np.sin(inc)],
            [0, np.sin(inc), np.cos(inc)]
        ])
        R_2_3 = np.array([
            [np.cos(AOP), -np.sin(AOP), 0],
            [np.sin(AOP), np.cos(AOP), 0],
            [0, 0, 1]
        ])

        R_0_3 = np.matmul(np.matmul(R_0_1,R_1_2),R_2_3)

        return R_0_3


    def kepler2cartesian(sma,ecc,inc, RAAN, AOP, TA):
        """This method converts Keplerian elements to Cartesian coordinates"""

        # Calculate the position in the orbital plane
        xyz = FrameTransformations.xyz_coordinates_in_plane(sma,ecc,TA)

        # Calculate the rotation matrix
        R = FrameTransformations.orbitPlane2Cartesian(inc, RAAN, AOP)

        # Rotate the position vector
        xyz = np.dot(R,xyz)
        return xyz

    def spherical2cartesian(r,theta,phi):
        """This method converts spherical coordinates to cartesian coordinates"""
        xyz = np.zeros(3)
        xyz[0] = r*np.cos(theta)*np.cos(phi)
        xyz[1] = r*np.sin(theta)*np.cos(phi)
        xyz[2] = r*np.sin(phi)
        return xyz

    def cartesian2spherical(xyz):
        """This method converts cartesian coordinates to spherical coordinates"""
        r = np.linalg.norm(xyz)
        phi = np.arcsin(xyz[2]/r)
        theta = np.arctan2(xyz[1],xyz[0])
        return r,theta,phi

    def eccentricAnomalyFromTrueAnomaly(ecc,ta):
        """This method calculates the eccentric anomaly from the true anomaly"""
        E = 2*np.arctan(np.tan(ta/2)/np.sqrt((1+ecc)/(1-ecc)))
        return E

    def trueAnomalyFromEccentricAnomaly(E,ecc):
        """This method calculates the true anomaly from the eccentric anomaly"""
        ta = 2*np.arctan(np.sqrt((1+ecc)/(1-ecc))*np.tan(E/2))
        return ta

    def meanAnomalyFromEccentricAnomaly(E,ecc):
        """This method calculates the mean anomaly from the eccentric anomaly"""
        M = E - ecc*np.sin(E)
        return M

    def eccentricAnomalyFromMeanAnomaly(M,ecc):
        """This method calculates the eccentric anomaly from the mean anomaly"""
        E = M
        while M-E+ecc*np.sin(E) > 1e-6:
            E = M + ecc*np.sin(E)
        return E

    def calculate_gmst(dateString):
        #"""
        #Calculate the Greenwich Mean Sidereal Time (GMST) for a given UTC date and time.
        #
        #Parameters:
        #date (datetime): The UTC date and time for which to calculate GMST.
        #
        #Returns:
        #float: The GMST in rad.
        #"""
        
        # Create an Astropy Time object
        t = Time(dateString, format='isot', scale='utc')

        # Get GMST in seconds
        gmst = t.sidereal_time('mean', 'greenwich')

        gmst_deg = gmst.deg  # Convert to degrees

        gmst_rad = gmst_deg*np.pi/180

        return gmst_rad

    def j2000_to_ecef(date_ref, xyz_j2000, elapsedTime):
    
        GMST_ref = FrameTransformations.calculate_gmst(date_ref)
    
        GMST = GMST_ref + EARTH_W*elapsedTime
    
        R_ecef_j2000 = np.array([
                        [np.cos(GMST), np.sin(GMST), 0],
                        [-np.sin(GMST), np.cos(GMST), 0],
                        [0,0,1]])
    
        xyz_ecef = np.dot(R_ecef_j2000, xyz_j2000)

        return xyz_ecef


if __name__ == "__main__":

    ecc = 0.00004650
    M0 = 117.052*np.pi/180
    E0 = FrameTransformations.eccentricAnomalyFromMeanAnomaly(M0, ecc)
    ta0 = FrameTransformations.trueAnomalyFromEccentricAnomaly(E0, ecc) 

    ta1 = np.arctan2(np.sin(E0)*np.sqrt(1-ecc**2),(np.cos(E0)-ecc))

    E1 = FrameTransformations.eccentricAnomalyFromTrueAnomaly(ecc, ta0)
    M1 = FrameTransformations.meanAnomalyFromEccentricAnomaly(E1,ecc)

    print(ta0*180/np.pi)
    print(ta1*180/np.pi)
  
    print(M1*180/np.pi)
    
    E2 = FrameTransformations.eccentricAnomalyFromTrueAnomaly(ecc, M0)
    M2 = FrameTransformations.meanAnomalyFromEccentricAnomaly(E2,ecc)

    print(M2*180/np.pi)