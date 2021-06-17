import xarray as xr
import numpy as np

import MDAnalysis as mda

from napalib.system.universe import NapAUniverse
from tqdm import tqdm

import os

def select_from_string(ag, string):
    atom, resid = string.split("_")
    atoms = ag.select_atoms(f"resid {resid} and name {atom}")
    assert len(atoms) == 1
    return atoms[0]

def collect(TOP, TRR, filename, dt=None, mutant=False):
    
    if os.path.exists(filename):
        print(f"{filename} already exists, please delete and run again")
        return -1
        
    u = NapAUniverse(TOP, mutant=True)
    u.load_new(TRR)
    
    N_frames = len(u.trajectory)
    
    if dt:
        times = np.arange(N_frames) * dt
    else:
        times = np.arange(N_frames)*u.trajectory.dt

    sodium_ag = u.select_atoms("resname SOD")
    
    A = u.atoms[u.atoms.segids == 'A']
    B = u.atoms[u.atoms.segids == 'B']

    atoms = [f"{ATOM}_{RESID}" for ATOM in ["OD1", "OD2"]
                              for RESID in [156, 157]]

    oxygens_a = mda.AtomGroup([select_from_string(A, a) for a in atoms])
    oxygens_b = mda.AtomGroup([select_from_string(B, a) for a in atoms])

    protomers = ['A', 'B']
    
    indices = np.zeros((N_frames, 2, 4), dtype=np.int32)
    distances = np.zeros((N_frames, 2, 4), dtype=np.float32)
    sodium = np.zeros((N_frames, len(sodium_ag), 2, 4), dtype=np.float32)
    
    for frame, ts in tqdm(enumerate(u.trajectory), total=N_frames):
        for pro, p in enumerate([oxygens_a, oxygens_b]):
            for i in range(4):
                reference_atom = p[i]
                vectors = reference_atom.position - sodium_ag.positions
                mdistances = np.linalg.norm(vectors, axis=1)
                smallest_distance_index = mdistances.argmin()
                distances[frame, pro, i] = mdistances[smallest_distance_index]
                indices[frame, pro, i] = sodium_ag[smallest_distance_index].index
                sodium[frame, :, pro, i] = mdistances
    
    data = xr.Dataset({
                        'distance': (('time', 'protomer', 'atom'), distances),
                        'sodium_index': (('time', 'protomer', 'atom'), indices),
                        'sodium_distances': (('time', 'sod_idx', 'protomer', 'atom'), sodium),
                        },
                      {'time': times,
                       'protomer': protomers,
                       'atom': atoms,
                       'sod_idx': sodium_ag.indices,
                       })
                
    data.to_netcdf(filename)
    return data
