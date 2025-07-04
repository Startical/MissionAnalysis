import numpy as np
import matplotlib.pyplot as plt


def plot_timeseries(timeseries, ylabel = '', units = '', title = '', labels = [], SF = 1):

    fig,ax = plt.subplots(1,1)

    Ndim = len(timeseries.data[0])

    if labels == []:
        labels = [f'Dim {i}' for i in range(Ndim)]

    for i in range(Ndim):
        data = [timeseries.data[j][i]*SF for j,t in enumerate(timeseries.time)]
        ax.plot(timeseries.time,data, label = labels[i])

    if units == '':
        units = timeseries.units

    ax.set_ylabel(f'{ylabel} {units}')
    ax.grid()

    ax.set_xlabel('Time (s)')
    ax.legend()

    if title == '':
        title = timeseries.id
    
    ax.set_title(f'{title}\n{timeseries.refTime}')

    fig.name = f'ts_{title}'

    return fig, ax


def plot_timeseries_vs_timeseries(timeseries, ylabel = '', units = '', title = '', labels = [], SF = 1):

    if not(isinstance(timeseries, list)):
        timeseries = [timeseries]

    Ndim = len(timeseries[0].data[0])

    if labels == []:
        labels = [timeseries[k].id for k in range(len(timeseries))]

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
            data = [ts.data[j][i]*SF for j,t in enumerate(time)]
            ax[i].plot(time,data, label = labels[k])

    if units == '':
        units = timeseries[0].units

    for i in range(Ndim):
        ax[i].set_ylabel(f'{ylabel} {i} {units}')
        ax[i].grid()

    ax[Ndim-1].set_xlabel('Time (s)')
    ax[Ndim-1].legend()

    if title == '':
        title = timeseries[0].id
    
    ax[0].set_title(f'{title}\n{timeseries[0].refTime}')

    fig.name = f'ts_vs_{title}'

    return fig, ax



