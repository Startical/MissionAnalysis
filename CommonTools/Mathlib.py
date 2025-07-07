import numpy as np
from scipy.spatial.transform import Rotation as R


def rotation_matrix_to_quat(R__A_B):
    """Convert a rotation matrix to a quaternion."""
    r = R.from_matrix(R__A_B)
    q__A_B = r.as_quat()  # returns (x, y, z, w) format)
    return q__A_B

def quat_mult(q__A_B, q__B_C):

    r__A_C = R.from_quat(q__A_B) * R.from_quat(q__B_C)
    q__A_C = r__A_C.as_quat()

    return q__A_C

def quat_conj(q__A_B):

    q__B_A = q__A_B

    q__B_A[0] = -q__B_A[0]
    q__B_A[1] = -q__B_A[1]
    q__B_A[2] = -q__B_A[2]

    return q__B_A

def quat_to_EulerAngles(q__A_B):

    r = R.from_quat(q__A_B)

    euler_zyx = r.as_euler('zyx')
    
    euler_rpy = euler_zyx[[2,1,0]]
    
    return euler_rpy

def quat_from_angle_vector(angle, vector):

    q = np.zeros(4)

    q[0] = np.sin(angle/2)*vector[0]
    q[1] = np.sin(angle/2)*vector[1]
    q[2] = np.sin(angle/2)*vector[2]
    
    q[3] = np.cos(angle/2)

    return q

if __name__ == "__main__":

    q_roll = quat_from_angle_vector(-30*np.pi/180, [1,0,0])
    q_pitch = quat_from_angle_vector(60*np.pi/180, [0,1,0])
    q_yaw = quat_from_angle_vector(45*np.pi/180, [0,0,1])

    q_rpy = quat_mult(q_roll,quat_mult(q_pitch, q_yaw))

    euler_rpy = quat_to_EulerAngles(q_rpy)

    print(euler_rpy*180/np.pi)
