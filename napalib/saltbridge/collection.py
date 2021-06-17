import MDAnalysis as mda
from napalib.system.universe import NapAUniverse
import numpy as np
import xarray as xr
from tqdm import tqdm

from .toptools import get_charged_residues, collection_scheme

import os


def collect(topology, trajectory, datafile_prefix, mutant=False):
    u = NapAUniverse(topology, mutant=mutant)
    u.load_new(trajectory)

    charged_residues = get_charged_residues(u)
    scheme = collection_scheme(charged_residues)

    N_frames = len(u.trajectory)
    dt = u.trajectory.dt
    time = [i * dt for i in range(N_frames)]

    # since the two titration states exist, there cannot be the same number of
    # salt bridges formed. This means that we need to store the data in two
    # separate arrays

    residues2str = lambda rg: [str(i) for i in rg]

    pos_A = list(filter(lambda x: x.positive, charged_residues['A']))
    neg_A = list(filter(lambda x: x.negative, charged_residues['A']))
    N_pos_A = len(pos_A)
    N_neg_A = len(neg_A)
    data_A = np.empty((N_frames, N_pos_A, N_neg_A), dtype=np.float32)
    da_A = xr.DataArray(data=data_A, dims=["time", "pos", "neg"],
                        coords=[time, residues2str(pos_A), residues2str(neg_A)])

    pos_B = list(filter(lambda x: x.positive, charged_residues['B']))
    neg_B = list(filter(lambda x: x.negative, charged_residues['B']))
    N_pos_B = len(pos_B)
    N_neg_B = len(neg_B)
    data_B = np.empty((N_frames, N_pos_B, N_neg_B), dtype=np.float32)
    da_B = xr.DataArray(data=data_B, dims=["time", "pos", "neg"],
                        coords=[time, residues2str(pos_B), residues2str(neg_B)])

    for frame, ts in tqdm(enumerate(u.trajectory), total=N_frames):
        for j, pair in enumerate(scheme['A']):
            i = j // N_neg_A
            k = j % N_neg_A
            da_A[frame, i, k] = pair.calculate_separation()
        for j, pair in enumerate(scheme['B']):
            i = j // N_neg_B
            k = j % N_neg_B
            da_B[frame, i, k] = pair.calculate_separation()
    # TODO these should really be writing to distancedata/
    da_A.to_netcdf(datafile_prefix + "_A.nc")
    da_B.to_netcdf(datafile_prefix + "_B.nc")

    return da_A, da_B
