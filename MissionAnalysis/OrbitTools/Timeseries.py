import numpy as np

class Timeseries(object):
    def __init__(self):
        self.time = []  # Time vector
        self.data = np.empty((0, 0))  # 2D data array

    def append(self, t, d):
        """Append a new time and data point"""
        self.time.append(t)
        if self.data.size == 0:
            self.data = np.array(d).reshape(1, -1)  # Initialize data array
        else:
            self.data = np.vstack([self.data, d])  # Append new data row

    def get_time(self):
        """Get the time vector"""
        return np.array(self.time)

    def get_data(self):
        """Get the data array"""
        return self.data

    def get_data_at_time(self, t):
        """Get the data at a specific time"""
        if t in self.time:
            index = self.time.index(t)
            return self.data[index]
        else:
            raise ValueError("Time not found in the time vector")

if __name__ == "__main__":
    # Example usage
    ts = Timeseries()
    ts.append(0, [1, 2, 3])
    ts.append(1, [4, 5, 6])
    ts.append(2, [7, 8, 9])

    print("Time vector:", ts.get_time())
    print("Data array:\n", ts.get_data())
    print("Data at time 1:", ts.get_data_at_time(1))




