from napalib.system.universe import NapAUniverse
from MDAnalysis.analysis.dihedrals import Dihedral
from MDAnalysis.core.groups import AtomGroup
import xarray as xr
import numpy as np
import os.path

from napalib.dihedrals.ags import *

def extract(TOP, TRAJ, cdf_fname):
    if os.path.exists(cdf_fname):
        print("{0} already exists, please delete and run again".format(cdf_fname))
        return -1

    u = NapAUniverse(TOP)
    u.load_new(TRAJ)
    prot = u.select_atoms("protein and backbone")
    protA = u.select_atoms("segid A")
    protB = u.select_atoms("segid B")

    phis_A = phi_pairs(protA.residues)
    psis_A = psi_pairs(protA.residues)

    ags_phi_A = list(map(get_phi_ag, phis_A))
    ags_psi_A = list(map(get_psi_ag, psis_A))
    ags_chi_1_A = list(map(get_chi_ag, protA.residues))[1:-1]

    ags_chi_1_A_numbered = [(i,j) for i,j in enumerate(ags_chi_1_A)]
    ags_chi_1_A_filtered = list(filter(lambda x: x[1], ags_chi_1_A_numbered))

    R_phi_A = Dihedral(ags_phi_A).run()
    R_psi_A = Dihedral(ags_psi_A).run()
    R_chi_A_reduced = Dihedral(list(map(lambda x: x[1] , ags_chi_1_A_filtered))).run()

    R_chi_A = np.zeros_like(R_phi_A.angles)
    for i,j  in enumerate(ags_chi_1_A_filtered):
        index = j[0]
        R_chi_A[:,index] = R_chi_A_reduced.angles[:,i]

    del R_chi_A_reduced

    pro_A_angles = np.stack([R_phi_A.angles, R_psi_A.angles, R_chi_A])
    angle, time, residue = pro_A_angles.shape
    # pro_A_angles = pro_A_angles.reshape((time, residue, angle))

    phis_B = phi_pairs(protB.residues)
    psis_B = psi_pairs(protB.residues)

    ags_phi_B = list(map(get_phi_ag, phis_B))
    ags_psi_B = list(map(get_psi_ag, psis_B))
    ags_chi_1_B = list(map(get_chi_ag, protB.residues))[1:-1]

    ags_chi_1_B_numbered = [(i,j) for i,j in enumerate(ags_chi_1_B)]
    ags_chi_1_B_filtered = list(filter(lambda x: x[1], ags_chi_1_B_numbered))

    R_phi_B = Dihedral(ags_phi_B).run()
    R_psi_B = Dihedral(ags_psi_B).run()
    R_chi_B_reduced = Dihedral(list(map(lambda x: x[1] , ags_chi_1_B_filtered))).run()

    R_chi_B = np.zeros_like(R_phi_B.angles)
    for i,j  in enumerate(ags_chi_1_B_filtered):
        index = j[0]
        R_chi_B[:,index] = R_chi_B_reduced.angles[:,i]

    del R_chi_B_reduced

    pro_B_angles = np.stack([R_phi_B.angles, R_psi_B.angles, R_chi_B])
    angle, time, residue = pro_B_angles.shape
    # pro_B_angles = pro_B_angles.reshape((time, residue, angle))

    stacked = np.stack([pro_A_angles, pro_B_angles])
    # stacked = xr.DataArray(np.stack([pro_A_angles, pro_B_angles]), 
    #                                          dims=['protomer', 'time', 'residue', 'angle'])
    stacked = xr.DataArray(np.stack([pro_A_angles, pro_B_angles]), 
                           dims=['protomer', 'angle', 'time', 'residue'])
    stacked = stacked.assign_coords(time=[i.time for i in u.trajectory], 
                                    protomer=["A", "B"], 
                                    angle=["phi", "psi", "chi_1"],
                                    residue=[i[0].resid for i in phis_B])
    stacked.to_netcdf(cdf_fname)
    return 0
