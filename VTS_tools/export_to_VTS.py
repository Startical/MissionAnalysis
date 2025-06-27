
import numpy as np
import os

from datetime import datetime

from Model.Spacecraft import Spacecraft
from OrbitTools.PlotUtils import initialize_ground_track_plot
import matplotlib.pyplot as plt

MJD0 = 55276
MJD_ref = datetime.strptime('2010-03-21T00:00:00Z', '%Y-%m-%dT%H:%M:%SZ')
user_name = os.getenv('USERNAME')

def get_default_destination_path():

    destination_path = os.path.join(os.path.dirname(__file__),'tmp')

    if not os.path.exists(destination_path):
        os.makedirs(destination_path)

    return destination_path


def get_time_reference_MJD_format(timeRef):
    '''
    VTS uses Modified Julian Date (MJD).
    '''

    time = datetime.strptime(timeRef, '%Y-%m-%dT%H:%M:%SZ')

    time_offset = (time-MJD_ref).total_seconds()

    MJD_day = np.floor(time_offset/(24*3600))
    MJD_seconds = time_offset - MJD_day * 24 * 3600

    MJD_day += MJD0  # Adjust to the VTS MJD reference

    return MJD_day, MJD_seconds

def export_position_file(xyz_j2000_ts, destination_path="", object_name = "MySat"):

    ## File destination and name
    if destination_path == "":
        destination_path = get_default_destination_path()

    filename = os.path.join(destination_path, f'{object_name}_position.txt')

    current_time = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")

    ## Time reference
    MJD_day, MJD_seconds = get_time_reference_MJD_format(xyz_j2000_ts.refTime)

    ## Open file for writing
    with open(filename, 'w') as f:
        # Write header
        f.write("CIC_OEM_VERS = 2.0\n")
        f.write(f"CREATION_DATE = {current_time}\n")
        f.write(f"ORIGINATOR = {user_name}\n")
        f.write("\n")
        f.write("META_START\n")
        f.write("\n")
        f.write(f"COMMENT Reference Date = {xyz_j2000_ts.refTime}\n")
        f.write("\n")
        f.write(f"OBJECT_NAME = {object_name}\n")
        f.write(f"OBJECT_ID = {object_name}\n")
        f.write("\n")
        f.write("CENTER_NAME = EARTH\n")
        f.write("REF_FRAME   = EME2000\n")
        f.write("TIME_SYSTEM = UTC\n")
        f.write("\n")
        f.write("META_STOP\n")
        f.write("\n")


        # Write data
        for i in range(len(xyz_j2000_ts.time)):
            t = xyz_j2000_ts.time[i]
            xyz = xyz_j2000_ts.data[i]

            tt = MJD_seconds+t
            days = np.floor(tt/24/3600)
            secs = tt - days*24*3600

            f.write(f"{MJD_day+days:.0f} {secs:.0f} {xyz[0]:.3f} {xyz[1]:.3f} {xyz[2]:.3f}\n")

    print(f"Position file exported to {filename}")

    return filename


def export_attitude_file(q_j2000_sc_ts, destination_path="", object_name = "MySat"):

    ## File destination and name
    if destination_path == "":
        destination_path = get_default_destination_path()

    filename = os.path.join(destination_path, f'{object_name}_attitude.txt')

    current_time = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")

    ## Time reference
    MJD_day, MJD_seconds = get_time_reference_MJD_format(q_j2000_sc_ts.refTime)

    ## Open file for writing
    with open(filename, 'w') as f:
        # Write header
        f.write("CIC_AEM_VERS = 1.0\n")
        f.write(f"CREATION_DATE = {current_time}\n")
        f.write(f"ORIGINATOR = {user_name}\n")
        f.write("\n")
        f.write("META_START\n")
        f.write("\n")
        f.write(f"COMMENT Reference Date = {q_j2000_sc_ts.refTime}\n")
        f.write("\n")
        f.write(f"OBJECT_NAME = {object_name}\n")
        f.write(f"OBJECT_ID = {object_name}\n")
        f.write("\n")
        f.write("REF_FRAME_A = EME2000\n")
        f.write("REF_FRAME_B = SC_BODY_1\n")
        f.write("ATTITUDE_DIR = A2B\n")
        f.write("\n")
        f.write("TIME_SYSTEM = UTC\n")
        f.write("\n")
        f.write("ATTITUDE_TYPE = QUATERNION\n")
        f.write("\n")
        f.write("META_STOP\n")
        f.write("\n")


        # Write data
        for i in range(len(q_j2000_sc_ts.time)):
            t = q_j2000_sc_ts.time[i]
            q = q_j2000_sc_ts.data[i]

            tt = MJD_seconds+t
            days = np.floor(tt/24/3600)
            secs = tt - days*24*3600

            f.write(f"{MJD_day+days:.0f} {secs:.0f} {q[3]:.6f} {q[0]:.6f} {q[1]:.6f} {q[2]:.6f}\n")

    print(f"Attitude file exported to {filename}")

    return filename


def export_spacecraft_pos_att_to_VTS(spacecraft, destination_path=""):

    export_position_file(spacecraft.xyz_j2000, destination_path, spacecraft.id)
    export_attitude_file(spacecraft.q_j2000_sc, destination_path, spacecraft.id)


if __name__ == "__main__":

    spacecraft = Spacecraft("IOD-1", kepler_parameters=np.array([7000, 0.01, np.radians(98), 0, 0, 0]), refTime='2025-03-25T12:00:00Z')

    spacecraft.propagate(90*60, 30)
    spacecraft.set_default_attitude_pointing(q_lof_sc = [0,0,0,1])

    export_spacecraft_pos_att_to_VTS(spacecraft, destination_path=r"C:\Users\plnegro\Programas\VTS\Data\IOD-1\Data")

    initialize_ground_track_plot(spacecraft)

    plt.show()
