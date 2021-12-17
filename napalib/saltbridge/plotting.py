import matplotlib.pyplot as plt
from .clustering import cluster_AP
import xarray as xr


def compare_clustering(file1, file2, stride=None, stride1=1, stride2=1):

    if stride:
        stride1 = stride
        stride2 = stride

    time1 = xr.load_dataarray(file1).time.values[::stride1] * 1e-6
    time2 = xr.load_dataarray(file2).time.values[::stride2] * 1e-6

    c1 = cluster_AP(file1, stride=stride1)
    c2 = cluster_AP(file2, stride=stride2)

    fig, axs = plt.subplots(2)

    ax1, ax2 = axs

    ax1.plot(time1, c1.labels_, '.', markersize=1)
    ax2.plot(time2, c2.labels_, '.', markersize=1)

    return fig, axs
