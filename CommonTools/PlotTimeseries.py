import numpy as np
import matplotlib.pyplot as plt


def plot_timeseries(timeseries, ylabel = '', units = '', title = ''):

    fig,ax = plt.subplots(1,1)

    Ndim = len(timeseries.data[0])

    for i in range(Ndim):
        data = [timeseries.data[j][i] for j,t in enumerate(timeseries.time)]
        ax.plot(timeseries.time,data, label = f'{i}')

    if units == '':
        units = timeseries.units

    ax.set_ylabel(f'{ylabel} {units}')
    ax.grid()

    ax.set_xlabel('Time (s)')
    ax.legend()

    if title == '':
        title = timeseries.id
    
    ax.set_title(title)

    return fig, ax


def plot_timeseries_vs_timeseries(timeseries, ylabel = '', units = '', title = ''):

    if not(isinstance(timeseries, list)):
        timeseries = [timeseries]

    Ndim = len(timeseries[0].data[0])

    fig,ax = plt.subplots(Ndim,1)

    # Convert axs to a flat list
    if isinstance(ax, np.ndarray):
        ax = ax.flatten().tolist()
    else:
        ax = [ax] # wrap single Axes object in a list

    for k in range(len(timeseries)):
        ts = timeseries[k]
        time = ts.time
        for i in range(Ndim):
            data = [ts.data[j][i] for j,t in enumerate(time)]
            ax[i].plot(time,data, label = timeseries[k].id)

    if units == '':
        units = timeseries[0].units

    for i in range(Ndim):
        ax[i].set_ylabel(f'{ylabel} {i} {units}')
        ax[i].grid()

    ax[Ndim-1].set_xlabel('Time (s)')
    ax[Ndim-1].legend()

    if title != '':
        ax[0].set_title(title)

    return fig, ax



