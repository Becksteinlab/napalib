from napalib.system.universe import NapAUniverse
from napalib.system.traj import Trajectory
import numpy as np
import xarray as xr
from tqdm import tqdm
from pathlib import Path

from .toptools import get_charged_residues, collection_scheme


def sep(atom1, atom2):
    return np.linalg.norm(atom1.position - atom2.position)


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

    def residues2str(rg):
        return [str(i) for i in rg]

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

    distancedata_dir = Path.cwd() / "distancedata"

    A_file = distancedata_dir / (datafile_prefix + "_A.nc")
    B_file = distancedata_dir / (datafile_prefix + "_B.nc")

    da_A.to_netcdf(A_file)
    da_B.to_netcdf(B_file)

    return da_A, da_B


def collect_305_156(trajectory: Trajectory):
    """Collect the K305-D157 salt bridge for both protomers and write out to a numpy file.

    This is potentially more useful than the standard collect function since it is protonation state independent.

    """

    u = trajectory.universe()

    NZ_A, NZ_B = u.select_atoms("resid 305 and name NZ")
    OD1_A, OD2_A, OD1_B, OD2_B = u.select_atoms("resid 156 and name OD1 OD2")

    data = np.zeros((3, u.trajectory.n_frames), dtype=np.float32)  # time, A, B

    for i, ts in tqdm(enumerate(u.trajectory), total=u.trajectory.n_frames):
        data[0, i] = ts.time

        distance_A = min(sep(NZ_A, OD1_A), sep(NZ_A, OD2_A))
        distance_B = min(sep(NZ_B, OD1_B), sep(NZ_B, OD2_B))

        data[1, i] = distance_A
        data[2, i] = distance_B

    datadir = Path.cwd() / "k305-d156-data"
    datadir.mkdir(exist_ok=True, parents=True)

    filename = datadir / f"{trajectory.name()}.npy"
    np.save(filename, data)


def collect_305_126(trajectory: Trajectory):
    """While not a salt bridge, we are generally interested in what happens to K305 when the K305-D156 salt bridge is
    broken.
    """

    u = trajectory.universe()

    NZ_A, NZ_B = u.select_atoms("resid 305 and name NZ")
    OG_A, OG_B = u.select_atoms("resid 126 and name OG1")

    data = np.zeros((3, u.trajectory.n_frames), dtype=np.float32)

    for i, ts in tqdm(enumerate(u.trajectory), total=u.trajectory.n_frames):
        data[0, i] = ts.time

        distance_A = sep(NZ_A, OG_A)
        distance_B = sep(NZ_B, OG_B)

        data[1, i] = distance_A
        data[2, i] = distance_B

    datadir = Path.cwd() / "k305-t126-data"
    datadir.mkdir(exist_ok=True, parents=True)

    filename = datadir / f"{trajectory.name()}.npy"
    np.save(filename, data)
