import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import seaborn as sns

import numpy as np
from numkit.timeseries import regularized_function

from napalib.theme import state_colors as colors


def stdev_histogrammed_function(t, y, **kwargs):
    return regularized_function(t, y, np.std, **kwargs)


def plot(RMSD_cdf, stride=1, ylim=[1.5, 11]):
    da = xr.load_dataarray(RMSD_cdf)
    time = da.time.values[::stride] * 1e-6
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))

    ax1, ax2 = axs

    ax1.plot(time, da.sel(crystal="OF", alignment="d", selection="c", protomer="A").values[::stride], label="OF")
    ax1.plot(time, da.sel(crystal="OCC", alignment="d", selection="c", protomer="A").values[::stride], label="OCC")
    ax1.plot(time, da.sel(crystal="IF", alignment="d", selection="c", protomer="A").values[::stride], label="IF")
    ax1.set_xlabel(r"time ($\mu$s)")
    ax1.set_ylabel(r"CRMSD ($\AA$)")
    ax1.set_ylim(ylim)

    ax2.plot(time, da.sel(crystal="OF", alignment="d", selection="c", protomer="B").values[::stride], label="OF")
    ax2.plot(time, da.sel(crystal="OCC", alignment="d", selection="c", protomer="B").values[::stride], label="OCC")
    ax2.plot(time, da.sel(crystal="IF", alignment="d", selection="c", protomer="B").values[::stride], label="IF")
    ax2.set_xlabel(r"time ($\mu$s)")
    ax2.set_ylim(ylim)

    ax2.legend(loc='best')

    return fig, axs


def publication_CRMSD_plot(RMSD_cdf, protomer, band_frames=100, xlim=None):
    da = xr.load_dataarray(RMSD_cdf)
    data_OF = da.sel(protomer=protomer,
                     selection='c',
                     alignment='d',
                     crystal="OF")

    data_IF = da.sel(protomer=protomer,
                     selection='c',
                     alignment='d',
                     crystal="IF")

    times = data_OF.time.values * 1e-6
    RMSD_OF = data_OF.values
    RMSD_IF = data_IF.values

    mavg_OF = np.convolve(RMSD_OF, np.ones(band_frames), 'valid') / band_frames
    mavg_IF = np.convolve(RMSD_IF, np.ones(band_frames), 'valid') / band_frames

    fig = plt.figure(figsize=(5, 2))
    ax = fig.add_subplot()

    sns.despine(ax=ax, offset=5)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['bottom'].set_linewidth(0.5)

    ax.set_xlabel(r"time ($\mu$s)")

    if xlim is None:
        ax.set_xlim([min(times), max(times)])

    ax.set_ylabel(r"CRMSD ($\AA$)")

    ax.yaxis.set_minor_locator(MultipleLocator(1))
    ax.xaxis.set_minor_locator(MultipleLocator(0.5))

    ax.plot(times[:-band_frames + 1], RMSD_OF[:-band_frames + 1], alpha=0.25, color=colors["OF"])
    ax.plot(times[:-band_frames + 1], mavg_OF, color=colors["OF"], label="5BZ3")

    ax.plot(times[:-band_frames + 1], RMSD_IF[:-band_frames + 1], alpha=0.25, color=colors["IF"])
    ax.plot(times[:-band_frames + 1], mavg_IF, color=colors["IF"], label="5BZ2")

    ax.legend(loc='best')

    plt.tight_layout()

    return fig, ax


def summary_plot(da, selection, alignment, xlims=[0.5, 2.5]):
    time = da.time * 1e-6
    A_IF = da.sel(protomer='A', crystal='IF',
                  selection=selection,
                  alignment=alignment).values
    A_OCC = da.sel(protomer='A', crystal='OCC',
                   selection=selection,
                   alignment=alignment).values
    A_OF = da.sel(protomer='A', crystal='OF',
                  selection=selection,
                  alignment=alignment).values
    B_IF = da.sel(protomer='B', crystal='IF',
                  selection=selection,
                  alignment=alignment).values
    B_OCC = da.sel(protomer='B', crystal='OCC',
                   selection=selection,
                   alignment=alignment).values
    B_OF = da.sel(protomer='B', crystal='OF',
                  selection=selection,
                  alignment=alignment).values

    fig, axs = plt.subplots(nrows=6, ncols=2, sharey='row')
    fig.set_size_inches(12, 12)

    color_if = "#7fc97f"
    color_occ = "#beaed4"
    color_of = "#fdc086"

    axs[0, 0].plot(time, A_IF, color=color_if)
    axs[1, 0].hist(A_IF, bins=50, density=True, color=color_if)

    axs[2, 0].plot(time, A_OCC, color=color_occ)
    axs[3, 0].hist(A_OCC, bins=50, density=True, color=color_occ)

    axs[4, 0].plot(time, A_OF, color=color_of)
    axs[5, 0].hist(A_OF, bins=50, density=True, color=color_of)

    axs[0, 1].plot(time, B_IF, color=color_if)
    axs[1, 1].hist(B_IF, bins=50, density=True, color=color_if)

    axs[2, 1].plot(time, B_OCC, color=color_occ)
    axs[3, 1].hist(B_OCC, bins=50, density=True, color=color_occ)

    axs[4, 1].plot(time, B_OF, color=color_of)
    axs[5, 1].hist(B_OF, bins=50, density=True, color=color_of)

    axs[0, 0].set_title("Protomer A (S1)")
    axs[0, 1].set_title("Protomer B (S4)")

    axs[0, 0].set_ylabel(r"IF RMSD ($\AA$)")
    axs[1, 0].set_ylabel("density")

    axs[2, 0].set_ylabel(r"OCC RMSD ($\AA$)")
    axs[3, 0].set_ylabel("density")

    axs[4, 0].set_ylabel(r"OF RMSD ($\AA$)")
    axs[5, 0].set_ylabel("density")

    if xlims is not None:
        for i in [1, 3, 5]:
            axs[i, 0].set_xlim(xlims)
            axs[i, 1].set_xlim(xlims)

    fig.tight_layout()

    return fig, axs
