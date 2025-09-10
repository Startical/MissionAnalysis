from skyfield.api import EarthSatellite
from OrbitTools.FrameTransformations import FrameTransformations as Frames
from Model.Spacecraft import Spacecraft

import numpy as np

from CommonTools.Datelib import time_offset

import requests

from tletools import TLE


def load_tle_from_file(filename, refDate, sat_id):
    '''
    The file is an export from https://www.space-track.org/#/gp 
    with multiple TLEs of the same object, sorted chronologically
    '''

    with open(filename, newline='', encoding='utf-8'):
        lines = [line.strip() for line in open(filename) if line.strip()]

    TLEs = {"data":[], "epoch":[]}
    for i in range(int(np.floor(len(lines)/2))):
        TLEs["data"].append("  "+sat_id+"\n"+lines[2*i]+"\n"+lines[2*i+1])
        epoch, *_ = tle2kepler(TLEs["data"][-1])
        TLEs["epoch"].append(epoch)

    tle_string = TLEs["data"][0]  # Use the first TLE for initialization

    for i in range(len(TLEs["epoch"])):

        if time_offset(TLEs["epoch"][i], refDate) < 0:
            tle_string = TLEs["data"][i]
        else: 
            break

    return tle_string

def get_TLE_with_norad_id(norad_id):
    # URL for the API endpoint
    celestrak_endpoint_url = "https://celestrak.org/NORAD/elements/gp.php?CATNR="
    # Make the GET request
    response = requests.get(celestrak_endpoint_url+norad_id, verify=False)

    # Check if the request was successful
    if response.status_code != 200:
        
        raise(ValueError(f'''
        Error retrieving TLE data: {response.status_code}
        Check Norad ID {norad_id}
        '''))

    return response.text


def fix_tle_spacing(tle_lines):
    # Remove leading blank spaces
    tle_lines = [line.lstrip() for line in tle_lines]

    # Ensure the first line has correct spacing
    line1 = tle_lines[0]
    line1 = line1[:1] + line1[1:7].rjust(6) + line1[7:8] + line1[8:17].rjust(9) + line1[17:20] + line1[20:32].rjust(12) + line1[32:43].rjust(11) + line1[43:51].rjust(8) + line1[51:60].rjust(9) + line1[60:68].rjust(8) + line1[68:]

    # Ensure the second line has correct spacing
    line2 = tle_lines[1]
    line2 = line2[:1] + line2[1:7].rjust(6) + line2[7:16].rjust(9) + line2[16:25].rjust(9) + line2[25:33].rjust(8) + line2[33:42].rjust(9) + line2[42:51].rjust(9) + line2[51:60].rjust(9) + line2[60:68].rjust(8) + line2[68:]

    return [line1, line2]

def tle2kepler(tle_string):
    spacecraft_id = tle_string.strip().splitlines()[0][2:]
    tle_lines = tle_string.strip().splitlines()[1:]  # skip the name line
    formatted_tle_lines = fix_tle_spacing(tle_lines)

    sat = EarthSatellite(formatted_tle_lines[0], formatted_tle_lines[1], spacecraft_id)

    epoch = sat.epoch.utc_iso()

    kepler_elements_v = [Frames.semiMajorAxisFromMeanMotion(sat.model.no_kozai/60),
                        sat.model.ecco,
                        sat.model.inclo,
                        sat.model.nodeo,
                        sat.model.argpo,
                        Frames.trueAnomalyFromMeanAnomaly(sat.model.mo,sat.model.ecco)]
    
    tle = TLE.from_lines(spacecraft_id,*formatted_tle_lines)
    orbit = tle.to_orbit()
    
    kepler_elements = orbit.classical()
    epoch2 = tle.epoch.strftime('%Y-%m-%dT%H:%M:%SZ')
    
    kepler_elements_v2 = [kepler_elements[0].value,
                        kepler_elements[1].value,
                        np.deg2rad(kepler_elements[2].value),
                        np.deg2rad(kepler_elements[3].value),
                        np.deg2rad(kepler_elements[4].value),
                        np.deg2rad(kepler_elements[5].value)]

    return epoch, kepler_elements_v, sat.model.satnum

def initialize_spacecraft_from_TLE(spacecraft_id, tle_string):

    epoch, kepler_elements_v, norad_id = tle2kepler(tle_string)

    spacecraft = Spacecraft(spacecraft_id, kepler_elements_v, epoch)
    spacecraft.norad_id = norad_id

    return spacecraft

if __name__ == '__main__':

    tle_string = """
    OBJECT C
    1 63212U 25052C   25092.84615115  .00046643  00000-0  22484-2 0  9991
    2 63212  97.4349 347.8550 0001428 176.9790 183.1453 15.18556553  2964
    """


    spacecraft = initialize_spacecraft_from_TLE("IOD1", tle_string)
    spacecraft.propagate(3600,60)

    print(spacecraft.id)
