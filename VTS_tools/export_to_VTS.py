
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

    first_line = "CIC_OEM_VERS = 2.0"
    comment_section = f"COMMENT Reference Date = {xyz_j2000_ts.refTime}\n"
    custom_section1 = f'''CENTER_NAME = EARTH
REF_FRAME   = EME2000'''

    header = build_header(first_line=first_line, object_name=object_name, custom_section1=custom_section1, comment_section=comment_section)

    data_string = timeseries_to_string(xyz_j2000_ts, format = '.3f')

    filename = export_file(destination_path=destination_path, file_name=f'{object_name}_position.txt', header=header, data=data_string)

    return filename

def export_attitude_file(q_j2000_sc_ts, destination_path="", object_name = "MySat"):
        
    first_line = "CIC_AEM_VERS = 1.0"
    comment_section = f"COMMENT Reference Date = {q_j2000_sc_ts.refTime}\n"
    custom_section1 = f'''REF_FRAME_A = EME2000
REF_FRAME_B = SC_BODY_1
ATTITUDE_DIR = A2B'''

    custom_section2= f'ATTITUDE_TYPE = QUATERNION'''

    header = build_header(first_line=first_line, object_name=object_name, custom_section1=custom_section1, custom_section2 = custom_section2, comment_section=comment_section)

    data_string = timeseries_to_string(q_j2000_sc_ts, order = [3,0,1,2], format = '.6f')
    
    filename = export_file(destination_path=destination_path, file_name=f'{object_name}_attitude.txt', header=header, data=data_string)

    return filename

def export_timeseries_file(ts, destination_path="", object_name = "MySat", order = [0,1,2,3,4], format = '.6f'):
        
    first_line = "CIC_MEM_VERS = 1.0"
    comment_section = f"COMMENT Reference Date = {ts.refTime}\n"
    custom_section1 = f'''USER_DEFINED_PROTOCOL = NONE
USER_DEFINED_CONTENT = {ts.id}
USER_DEFINED_SIZE = {len(ts.data[0])}
USER_DEFINED_TYPE = REAL
USER_DEFINED_UNIT = {ts.units}'''

    header = build_header(first_line=first_line, object_name=object_name, custom_section1=custom_section1, comment_section=comment_section)

    data_string = timeseries_to_string(ts, order=order, format=format)
    
    filename = export_file(destination_path=destination_path, file_name=f'{object_name}_{ts.id}.txt', header=header, data=data_string)

    return filename

def build_header(first_line = 'CIC_MEM_VERS = 1.0', object_name = '', custom_section1 = '', custom_section2 = '', comment_section=''):

    header = first_line + '\n'

    current_time = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")

    header = header + f"CREATION_DATE = {current_time}\n"
    header = header + f"ORIGINATOR = {user_name}\n"
    header = header + "\n"
    header = header + "META_START\n"
    header = header + "\n"
    header = header + comment_section+ "\n"
    header = header + f"OBJECT_NAME = {object_name}\n"
    header = header + f"OBJECT_ID = {object_name}\n"
    header = header + "\n"
    header = header + custom_section1 + "\n"
    header = header + "TIME_SYSTEM = UTC\n"
    header = header + "\n"
    header = header + custom_section2 + "\n"
    header = header + "\n"
    header = header + "META_STOP\n"
    header = header + "\n"

    return header

def timeseries_to_string(ts, order = [0,1,2,3,4], format = '.6f'):

    ## Time reference
    MJD_day, MJD_seconds = get_time_reference_MJD_format(ts.refTime)

    # Write data
    ts_string = ''
    for i in range(len(ts.time)):
        t = ts.time[i]
        data = ts.data[i]

        tt = MJD_seconds+t
        days = np.floor(tt/24/3600)
        secs = tt - days*24*3600

        ts_line = f"{MJD_day+days:.0f} {secs:.0f}"
        
        for j in range(len(data)):
            ts_line = ts_line + f" {data[order[j]]:{format}}"
            
        ts_string = ts_string + ts_line + "\n"

    return ts_string

def export_file(destination_path="", file_name="my_file.txt", header="", data=""):
    ## File destination and name
    if destination_path == "":
        destination_path = get_default_destination_path()
    filename = os.path.join(destination_path, file_name)
    ## Open file for writing
    with open(filename, 'w') as f:
        # Write header
        f.write(header)
        
        # Write data
        f.write(data)
    print(f"File exported to {filename}")
    return filename




def export_spacecraft_pos_att_to_VTS(spacecraft, destination_path=""):

    export_position_file(spacecraft.xyz_j2000, destination_path, spacecraft.id)
    export_attitude_file(spacecraft.q_j2000_sc, destination_path, spacecraft.id)


if __name__ == "__main__":

    spacecraft = Spacecraft("IOD-1", kepler_parameters=np.array([7000, 0.01, np.radians(98), 0, 0, 0]), refTime='2025-03-25T12:00:00Z')

    spacecraft.propagate(90*60, 30)
    spacecraft.set_default_attitude_pointing(q_lof_sc = [0,0,0,1])

    export_spacecraft_pos_att_to_VTS(spacecraft, destination_path=r"C:\Users\plnegro\Programas\VTS\Data\IOD-1\Data")
    export_timeseries_file(spacecraft.ta, destination_path=r"C:\Users\plnegro\Programas\VTS\Data\IOD-1\Data", object_name=spacecraft.id, format='.3f')


    initialize_ground_track_plot(spacecraft)

    plt.show()
