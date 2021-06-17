import xarray as xr
import numpy as np

import MDAnalysis as mda
from MDAnalysis.analysis import rms
from MDAnalysis.analysis import align
import os

from napalib.system.universe import NapAUniverse


def get_CRMSD(ref_align_ag, mobile_align_ag, ref_select_ag, mobile_select_ag):
    ref0 = ref_align_ag.positions - ref_align_ag.center_of_mass()
    mob0 = mobile_align_ag.positions - mobile_align_ag.center_of_mass()
    R, rmsd = align.rotation_matrix(mob0, ref0)

    both = mobile_align_ag + mobile_select_ag

    both.translate(-mobile_align_ag.center_of_mass())
    both.rotate(R)
    both.translate(ref_align_ag.center_of_mass())

    rmsd = rms.rmsd(mobile_select_ag.positions, ref_select_ag.positions)
    return rmsd


def extract_RMSD(TOP, TRR, filename, mutant=False):
    if os.path.exists(filename):
        print("% already exists, please delete and run again" % filename)
        return -1

    from napalib.system.ref import IF_crystal_ref as uIF
    from napalib.system.ref import OF_crystal_ref as uOF
    from napalib.system.ref import OCC_ref as uOCC
    from napalib.system.ref import IF_ref as uEIF

    u = NapAUniverse(TOP, mutant=mutant)
    u.load_new(TRR)
    N_frames = len(u.trajectory)

    data = np.zeros((3, 3, 2, 4, N_frames), dtype=np.float32)

    dt = u.trajectory[1].time - u.trajectory[0].time
    time = [i * dt for i in range(len(u.trajectory))]

    da = xr.DataArray(data, dims=['alignment',
                                  'selection',
                                  'protomer',
                                  'crystal',
                                  'time'])

    da = da.assign_coords(alignment=['c', 'd', 'f'],
                          selection=['c', 'd', 'f'],
                          protomer=['A', 'B'],
                          crystal=["OF", "IF", "OCC", "EIF"],
                          time=time)

    ag_A = u.select_atoms("backbone and segid A")
    ag_B = u.select_atoms("backbone and segid B")

    ag_A_c = u.core_A.select_atoms("backbone")
    ag_A_d = u.dimer_A.select_atoms("backbone")
    ag_A_f = ag_A_c + ag_A_d + u.neither_A.select_atoms("backbone")

    ag_B_c = u.core_B.select_atoms("backbone")
    ag_B_d = u.dimer_B.select_atoms("backbone")
    ag_B_f = ag_B_c + ag_B_d + u.neither_B.select_atoms("backbone")

    assert len(ag_A) == len(ag_B)

    ag_of_xtal_c = uOF.core_A.select_atoms("backbone and not bynum 4921 2221")
    ag_of_xtal_d = uOF.dimer_A.select_atoms("backbone and not bynum 4921 2221")
    ag_of_xtal_f = ag_of_xtal_c + ag_of_xtal_d + uOF.neither_A.select_atoms("backbone and not bynum 4921 2221")

    ag_if_xtal_c = uIF.core_A.select_atoms("backbone")
    ag_if_xtal_d = uIF.dimer_A.select_atoms("backbone")
    ag_if_xtal_f = ag_if_xtal_c + ag_if_xtal_d + uIF.neither_A.select_atoms("backbone")

    ag_occ_c = uOCC.core_A.select_atoms("backbone")
    ag_occ_d = uOCC.dimer_A.select_atoms("backbone")
    ag_occ_f = ag_occ_c + ag_occ_d + uOCC.neither_A.select_atoms("backbone")

    ag_eif_c = uEIF.core_A.select_atoms("backbone")
    ag_eif_d = uEIF.dimer_A.select_atoms("backbone")
    ag_eif_f = ag_eif_c + ag_eif_d + uEIF.neither_A.select_atoms("backbone")

    msas = {}
    msas["A"] = {'c': ag_A_c, 'd': ag_A_d, 'f': ag_A_f}
    msas["B"] = {'c': ag_B_c, 'd': ag_B_d, 'f': ag_B_f}

    rsas = {}
    rsas["OF"] = {'c': ag_of_xtal_c, 'd': ag_of_xtal_d, 'f': ag_of_xtal_f}
    rsas["IF"] = {'c': ag_if_xtal_c, 'd': ag_if_xtal_d, 'f': ag_if_xtal_f}
    rsas["OCC"] = {'c': ag_occ_c, 'd': ag_occ_d, 'f': ag_occ_f}
    rsas["EIF"] = {'c': ag_eif_c, 'd': ag_eif_d, 'f': ag_eif_f}

    maas = msas
    raas = rsas

    protein_AG = u.select_atoms("protein and name CA")

    print(f"Collecting {TOP}")

    for frame, ts in enumerate(u.trajectory):
        for alignment in da.coords["alignment"].data:
            for selection in da.coords["selection"].data:
                for protomer in da.coords["protomer"].data:
                    for crystal in da.coords["crystal"].data:
                        mobile_select_ag = msas[protomer][selection]
                        ref_select_ag = rsas[crystal][selection]
                        mobile_align_ag = maas[protomer][alignment]
                        ref_align_ag = raas[crystal][alignment]
                        assert len(ref_align_ag) == len(mobile_align_ag)

                        ref0 = ref_align_ag.positions - ref_align_ag.center_of_mass()
                        mob0 = mobile_align_ag.positions - mobile_align_ag.center_of_mass()
                        R, rmsd = align.rotation_matrix(mob0, ref0)

                        prot = {"A": ag_A, "B": ag_B}

                        prot[protomer].translate(-mobile_align_ag.center_of_mass())
                        prot[protomer].rotate(R)
                        prot[protomer].translate(ref_align_ag.center_of_mass())

                        rmsd = rms.rmsd(mobile_select_ag.positions, ref_select_ag.positions)

                        da.sel(alignment=alignment,
                               selection=selection,
                               protomer=protomer,
                               crystal=crystal)[frame] = rmsd
    da.to_netcdf(filename)
    return 0
