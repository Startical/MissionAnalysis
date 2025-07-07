import numpy as np


class Timeseries(object):

    id = ""  # Identifier for the timeseries
    refTime = 'YYYY-MM-DDTHH:MM:SSZ'
    units = "[n/a]"  # Units of the timeseries data

    def __init__(self, ts_id = "", refTime = 'YYYY-MM-DDTHH:MM:SSZ'):
        self.id = ts_id  # Identifier for the timeseries
        self.time = []  # Time vector
        self.data = np.empty((0, 0))  # 2D data array
        self.refTime = refTime  # Reference time for the timeseries

    def append(self, t, d):
        """Append a new time and data point"""
        self.time.append(t)
        if self.data.size == 0:
            self.data = np.array(d).reshape(1, -1)  # Initialize data array
        else:
            self.data = np.vstack([self.data, d])  # Append new data row

    def cleanup_timeseries(self):
        """Sort in increasing vaues of time and remove duplicated times"""
        # sort in increasing values of time
        idx = [i for i, _ in sorted(enumerate(self.time), key=lambda x: x[1])]

        self.time = [self.time[i] for i in idx]
        self.data = self.data[idx]

        # remove duplicated items
        
        # Get indexes of first occurrences of unique values
        unique_indexes = []
        seen = set()
        for i, val in enumerate(self.time):
            if val not in seen:
                seen.add(val)
                unique_indexes.append(i)

        self.time = [self.time[i] for i in unique_indexes]
        self.data = self.data[unique_indexes]

    def get_time(self):
        """Get the time vector"""
        return np.array(self.time)

    def get_data(self):
        """Get the data array"""
        return self.data

    def get_index(self, t):
        """Get the index of a specific time"""
        if t in self.time:
            index = self.time.index(t)
            return index
        else:
            raise ValueError("Time not found in the time vector")


    def get_data_at_time(self, t):
        """Get the data at a specific time"""
        if t in self.time:
            index = self.time.index(t)
            return self.data[index]
        elif t < self.time[0] or t > self.time[-1]:
            raise ValueError("Time is outside the range of the time vector")
        else:
            # Find the two nearest time points
            for i in range(len(self.time) - 1):
                t0, t1 = self.time[i], self.time[i + 1]
                if t0 < t < t1:
                    d0, d1 = self.data[i], self.data[i + 1]
                    # Linear interpolation
                    weight = (t - t0) / (t1 - t0)
                    return d0 + weight * (d1 - d0)
            raise ValueError("Interpolation failed: time not bracketed properly")


    def get_data_at_index(self, idx):
        """Get the data at a specific index"""
        if idx < len(self.time):
            return self.time[idx], self.data[idx]
        else:
            raise ValueError("Index out of bounds")

if __name__ == "__main__":
    # Example usage
    ts = Timeseries()
    ts.append(0, [1, 2, 3])
    ts.append(1, [4, 5, 6])
    ts.append(2, [7, 8, 9])

    print("Time vector:", ts.get_time())
    print("Data array:\n", ts.get_data())
    print("Data at time 1:", ts.get_data_at_time(1))
    print("Data at time (last):", ts.get_data_at_index(-1))



