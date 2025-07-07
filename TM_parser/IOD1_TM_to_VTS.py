from TM_parser.IOD1_TM_parser import *
from VTS_tools.export_to_VTS import *

import CommonTools.Datelib as Dates
import CommonTools.Mathlib as Math
from CommonTools.Timeseries import Timeseries

import os
from Model.SatelliteUtils import initialize_spacecraft_from_TLE
from Model.SatelliteUtils import tle2kepler
from CommonTools.PlotTimeseries import *


def export_TM_and_plots(TM_file_name, tles, save_figs = False, results_folder = ''):

    if not os.path.exists(results_folder):
        os.makedirs(results_folder)

    ts = read_timeseries_from_file(TM_file_name)


    figs = ['','','','','','','']
    
    figs[0],_ = plot_timeseries(ts['ctrl_refq'], labels = ['qx', 'qy', 'qz', 'q0'])
    figs[1],_ = plot_timeseries(ts['ctrl_errq'], labels = ['qx', 'qy', 'qz', 'q0'])

    tle_string = tles["data"][0]  # Use the first TLE for initialization

    for i in range(len(tles["epoch"])):

        if Dates.time_offset(tles["epoch"][i], ts['ctrl_refq'].refTime) < 0:
            tle_string = tles["data"][i]
        else: 
            break

    spacecraft = initialize_spacecraft_from_TLE("IOD-1", tle_string)
    spacecraft.propagate_from_start_date(ts['ctrl_refq'].refTime,3600,1)

    spacecraft.set_default_attitude_pointing(q_lof_sc = [0,0,0,1])

    q_j2000_sc_ref_ts = spacecraft.q_j2000_sc
    q_j2000_sc_ref_ts.id = 'q_j2000_sc_ref'
    q_j2000_lof_ts = spacecraft.q_j2000_lof

    q_j2000_sc_ts = Timeseries('q_j2000_sc',ts['ctrl_refq'].refTime)
    q_lof_sc_ts = Timeseries('q_lof_sc',ts['ctrl_refq'].refTime)
    eulerAngles_ts = Timeseries('eulerAngles',ts['ctrl_refq'].refTime)

    for i,t in enumerate(ts['ctrl_refq'].time):

        qRef = ts['ctrl_refq'].get_data_at_time(t)

        t1 = Dates.time_offset(ts['ctrl_errq'].refTime, ts['ctrl_refq'].refTime) + t
        qErr = ts['ctrl_errq'].get_data_at_time(t1)

        q_j2000_sc = Math.quat_mult(qRef, qErr)
        q_j2000_sc_ts.append(t,q_j2000_sc)

        t2 = Dates.time_offset(q_j2000_lof_ts.refTime, ts['ctrl_refq'].refTime) + t
        q_j2000_lof = q_j2000_lof_ts.get_data_at_time(t2)

        q_lof_sc = Math.quat_mult(Math.quat_conj(q_j2000_lof), q_j2000_sc)
        q_lof_sc_ts.append(t,q_lof_sc)
        
        eulerAngles = Math.quat_to_EulerAngles(q_lof_sc)
        eulerAngles_ts.append(t,eulerAngles)

    figs[2],_ = plot_timeseries_vs_timeseries([ts['ctrl_refq'],q_j2000_sc_ref_ts])
    figs[3],_ = plot_timeseries_vs_timeseries([q_j2000_sc_ts,q_j2000_sc_ref_ts])
    figs[4],_ = plot_timeseries(q_lof_sc_ts, labels = ['qx', 'qy', 'qz', 'q0'])
    figs[5],_ = plot_timeseries(eulerAngles_ts, ylabel = 'Angle', units = '(deg)', labels = ['roll', 'pitch', 'yaw'], SF = 180/np.pi)


    export_position_file(spacecraft.xyz_j2000, destination_path=results_folder, object_name=spacecraft.id)
    export_attitude_file(q_j2000_sc_ts, destination_path=results_folder, object_name=spacecraft.id)
    export_timeseries_file(eulerAngles_ts, destination_path=results_folder, object_name=spacecraft.id, format='.3f')

    figs[6],_ = initialize_ground_track_plot(spacecraft)

    if save_figs:
        for fig in figs:
            fig_name = os.path.join(results_folder, f"{fig.name}.png")
            fig.savefig(fig_name)
            print(f"Figure saved: {fig_name}")

        plt.close('all')
    else:
        plt.show()

def load_tles(filename):

    with open(filename, newline='', encoding='utf-8'):
        lines = [line.strip() for line in open(filename) if line.strip()]

    TLEs = {"data":[], "epoch":[]}
    for i in range(int(np.floor(len(lines)/2))):
        TLEs["data"].append("  IOD1\n"+lines[2*i]+"\n"+lines[2*i+1])
        epoch, *_ = tle2kepler(TLEs["data"][-1])
        TLEs["epoch"].append(epoch)

    return TLEs



if __name__ == "__main__":

    
    dir_name = os.path.dirname(__file__)
    results_folder = dir_name #r"C:\Users\plnegro\Programas\VTS\Data\IOD-1\Data"
    TM_folder = "TM_data"

    tests = ["16Jun_day","16Jun_night","20Jun_day","23Jun_night","24Jun_night","26Jun_day","26Jun_night","27Jun_day","27Jun_night","28Jun_day","30Jun_night" ]
    

    tles = load_tles(os.path.join(dir_name,"TLEs.txt"))

    save_figs = True

    for test in tests:
        TM_file_name = os.path.join(dir_name,TM_folder, f"test_quats_{test}.csv")

        results_folder = os.path.join(dir_name,f"test_{test}")

        export_TM_and_plots(TM_file_name, tles, save_figs = save_figs ,results_folder = results_folder)


