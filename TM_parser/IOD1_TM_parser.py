import numpy as np

import csv
import re

from CommonTools.Timeseries import Timeseries
import CommonTools.Datelib as Dates


def read_timeseries_from_file(filename, header_lines=1):

    timeseries = {}

    with open(filename, newline='', encoding='utf-8') as csvfile:

        reader = csv.reader(csvfile)

        for index, row in enumerate(reader):
            if index < header_lines:
                continue ## skip the header rows

            ## Parse
            dateRaw = row[5]
            date = re.sub(r'\.(\d{3})\d*(?=Z)', r'.\1', dateRaw)

            parameterId = row[4]

            data = []
            value_ok = True

            try:
                for i in range(len(row)-6):

                    value = row[6+i]
                
                    if data == 'NaN':
                        break
                    else:
                        data.append(float(value))
            except:
                value_ok = False

            if value_ok:

                ## Assign to timeseries
                if parameterId in timeseries:
                    ts = timeseries[parameterId]
                else:
                    ts = Timeseries(parameterId,date)
                    timeseries[parameterId] = ts

                t =  Dates.time_offset(date,ts.refTime)

                ts.append(t, data)

    for ts_id in timeseries.keys():
        timeseries[ts_id].cleanup_timeseries()

    return timeseries

if __name__ == "__main__":

    
    import os
    dir_name = os.path.dirname(__file__)
    file_name = os.path.join(os.path.dirname(__file__),"test_quats_27Jun.csv")

    ts = read_timeseries_from_file(file_name)

    from CommonTools.PlotTimeseries import *

    plot_timeseries(ts['ctrl_refq'])
    plot_timeseries(ts['ctrl_errq'])

    #plot_timeseries_vs_timeseries([ts['bias_ok'], ts['bias_ok']])

    plt.show()





