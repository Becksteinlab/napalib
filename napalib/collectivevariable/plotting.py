import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import xarray as xr
import math

import seaborn as sns

reference = {"OF" : {"dz" : 6.93, "dphi" : 0},
             "IF" : {"dz" : -1.65, "dphi" : 11.92}
             }

def publication_timeseries_plot(ncdf, protomer, band_frames=100, with_ref=True):
    data = xr.load_dataset(ncdf)
    
    z_data = data['dz'].sel(protomer=protomer).values
    phi_data = data['dphi'].sel(protomer=protomer).values
    time = data.time.values * 1e-6
    
    fig = plt.figure()
    gs = fig.add_gridspec(2, 2, hspace=0.1, width_ratios=[3,1])
    (phi, hist_phi), (z, hist_z) = gs.subplots(sharex=False, sharey=False)

    sns.despine(ax=phi, offset=10)
    sns.despine(ax=hist_phi, bottom=True)
    sns.despine(ax=z, offset=10)
    sns.despine(ax=hist_z, bottom=True)

    # phi axis tweaks
    degree_formatter = plt.matplotlib.ticker.StrMethodFormatter(r"{x:g}$\degree$")
    phi.yaxis.set_major_formatter(degree_formatter)
    phi.axes.xaxis.set_visible(False)
    phi.set_ylabel("$\Delta \phi$")
    plt.setp(phi.get_xticklabels(), visible=False)
    phi.yaxis.set_minor_locator(MultipleLocator(2))
    phi.set_xlim([0, math.ceil(max(time))])
    phi.spines['bottom'].set_visible(False)
    phi.spines['left'].set_linewidth(0.5)
    hist_phi.spines['left'].set_linewidth(0.5)
    hist_phi.axes.xaxis.set_visible(False)
    hist_phi.yaxis.set_major_formatter(degree_formatter)

    mavg_phi = np.convolve(phi_data, np.ones(band_frames), 'valid') / band_frames
    phi.plot(time, phi_data, alpha=0.25, color='k')
    phi.plot(time[:-band_frames+1], mavg_phi, color='k')
    phi.axhline(reference['OF']['dphi'])
    phi.axhline(reference['IF']['dphi'])
    
    hist_phi.hist(phi_data , histtype='step', bins=40, density=True, orientation='horizontal', color='k')

    # z axis tweaks
    z.set_ylabel("$\Delta$Z ($\AA$)")
    z.set_xlabel("time ($\mu$s)")
    z.yaxis.set_minor_locator(MultipleLocator(1))
    hist_z.axes.xaxis.set_visible(False)
    z.spines['left'].set_linewidth(0.5)
    z.spines['bottom'].set_linewidth(0.5)
    hist_z.spines['left'].set_linewidth(0.5)

    mavg_z = np.convolve(z_data, np.ones(band_frames), 'valid') / band_frames
    z.plot(time, z_data, alpha=0.25, color='k')
    z.plot(time[:-band_frames+1], mavg_z, color='k')
    z.axhline(reference['OF']['dz'])
    z.axhline(reference['IF']['dz'])

    hist_z.hist(z_data, histtype='step', bins=40, density=True, orientation='horizontal', color='k')

    fig.subplots_adjust(bottom=0.150)

    return fig, (phi, z)

def publication_hist_plot(ncdf, protomer, ax=None, with_ref=True, **kwargs):
    
    kwargs.setdefault("color", "crimson")
    data = xr.load_dataset(ncdf)
    
    if ax is None:
        fig = plt.figure(figsize=(4,3))
        ax = fig.add_subplot()
    ax = plot_hist(data.sel(protomer=protomer), ax, **kwargs)
    ax.plot([reference["OF"]["dphi"]],[reference["OF"]["dz"]], 's', color='k', label="OF")
    ax.plot([reference["IF"]["dphi"]],[reference["IF"]["dz"]], 'v', color='k', label="IF")
    ax.legend(loc='best')
    
    ax.set_xlabel("$\Delta \phi$")
    ax.set_ylabel("$\Delta$Z ($\AA$)")
    
    ax.yaxis.set_minor_locator(MultipleLocator(1))
    ax.xaxis.set_minor_locator(MultipleLocator(2))
    degree_formatter = plt.matplotlib.ticker.StrMethodFormatter(r"{x:g}$\degree$")
    ax.xaxis.set_major_formatter(degree_formatter)
    sns.despine(ax=ax)
    
    fig = ax.figure
    fig.tight_layout()
    return fig, ax


def plot_hist(dataset, ax, color, bins=50, alpha=0.2,
              frac_levels=[0.98, 0.5, 0.25],
              phi_range=[-5, 45], z_range=[-12, 12], filled=True):
    z = dataset['dz'].values
    phi = dataset['dphi'].values
    h, x, y = np.histogram2d(phi, z, bins=bins, range=[phi_range, z_range])
    h_norm = h/np.sum(h)
    h_vals = [i for i in np.sort(h_norm, axis=None) if i != 0]
    frac_counts = [sum(h_vals[i:]) for i in range(len(h_vals))]
    for fi, flevel in enumerate(frac_levels):
        vlevel = h_vals[
                  [i for i, val in enumerate(frac_counts) if val <= flevel][0]]
        ax.contour(h_norm.T, levels=[vlevel-1e-5, 1], colors=color,
                   alpha=alpha, extent=[*phi_range, *z_range], linewidths=1)
        if filled:
            ax.contourf(h_norm.T, levels=[vlevel-1e-5, 1], colors=color,
                        alpha=alpha, extent=[*phi_range, *z_range])
    return ax


def quick_hist(dataset, color='crimson', **kwargs):
    fig = plt.figure(num=None, figsize=(6.5, 3))
    fig.subplots_adjust(bottom=0.20, wspace=0)
    A, B = fig.subplots(1, 2)

    A.annotate('A', (0.925, 0.05), xycoords='axes fraction', style="italic")
    B.annotate('B', (0.925, 0.05), xycoords='axes fraction', style="italic")

    A.set_xlabel(r"$\Delta \phi$ ($\degree$)")
    A.set_ylabel(r"$\Delta Z$ ($\AA$)")
    B.set_xlabel(r"$\Delta \phi$ ($\degree$)")

    A.yaxis.set_minor_locator(MultipleLocator(1))
    A.xaxis.set_minor_locator(MultipleLocator(2))
    B.yaxis.set_minor_locator(MultipleLocator(1))
    B.xaxis.set_minor_locator(MultipleLocator(2))

    B.set_yticklabels([])

    plot_hist(dataset.sel(protomer='A'), A, color, **kwargs)
    plot_hist(dataset.sel(protomer='B'), B, color, **kwargs)
    return fig, (A, B)
