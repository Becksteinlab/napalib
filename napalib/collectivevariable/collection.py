import numpy as np
from napalib.system.universe import NapAUniverse
from MDAnalysis.analysis import align
from MDAnalysis.analysis.leaflet import LeafletFinder
import xarray as xr
import os

from tqdm import tqdm


def get_z(core_ag, dimer_ag):
    """Returns the difference between the core and dimerisation domain
    centers of mass.

    Parameters:
    core_ag (AtomGroup): Atom group containing the core CA atoms
    dimer_ag (AtomGroup): Atom group containing the dimer CA atoms
    """
    return core_ag.center_of_mass()[2] - dimer_ag.center_of_mass()[2]


def get_z_memb(domain_ag, memb_ag):
    """Returns the difference between the dimerisation domain and membrane
    centers of mass.

    Parameters:
    domain_ag (AtomGroup): Atom group containing the domain CA atoms
    memb_ag (AtomGroup): Atom group containing the membrane atoms
    """
    return domain_ag.center_of_mass()[2] - memb_ag.center_of_mass()[2]


def get_dz(dimer_bundle_va, core_bundle_va):
    """Return the difference between the two virtual atoms of the dimer
    bundle and the core bundle.

    Parameters:
    dimer_bundle_va (numpy array): Position of the virtual atom defining
                                   the dimer bundle
    core_bundle_va (numpy array): Position of the virtual atom defining
                                  the core bundle
    """
    return core_bundle_va[2] - dimer_bundle_va[2]


def get_theta(top_core_va, bot_core_va, bot_dimer_va):
    """Returns the delta-theta angle.

    Parameters:
    top_core_va (numpy array): Position of the virtual atom defining the
                               top of the core domain.
    bot_core_va (numpy array): Position of the virtual atom defining the
                               bottom of the core domain.
    bot_dimer_va (numpy array): Position of the virtual atom defining the
                                bottom of the dimerization domain.
    """

    a_1 = top_core_va
    a_2 = bot_core_va
    a_3 = bot_dimer_va

    a_12 = a_1 - a_2
    a_32 = a_3 - a_2

    return np.rad2deg(np.arccos((a_12 / np.linalg.norm(a_12)).dot(a_32 / np.linalg.norm(a_32))))


def get_phi(ref_dimer_ag, ref_core_ag, dimer_ag, core_ag):
    """Returns the delta-phi angle. Note that this operation changes positions
    of the atoms involved--state is altered!

    Parameters:
    ref_dimer_ag (AtomGroup): Atom group containing the reference structure's
                              dimerisation domain CA atoms.
    ref_core_ag (AtomGroup): Atom group containing the reference structure's
                             core domain CA atoms.
    dimer_ag (AtomGroup): Atom group containing the mobile structure's
                          dimerisation domain CA atoms.
    core_ag (AtomGroup): Atom group containing the mobile structure's core
                         domain CA atoms.
    """
    ref0 = ref_dimer_ag.positions - ref_dimer_ag.center_of_mass()
    mob0 = dimer_ag.positions - dimer_ag.center_of_mass()
    R, rsmd = align.rotation_matrix(mob0, ref0)
    both = dimer_ag + core_ag
    both.translate(-dimer_ag.center_of_mass())
    both.rotate(R)
    both.translate(ref_dimer_ag.center_of_mass())

    R, rmsd = align.rotation_matrix(core_ag.positions - core_ag.center_of_mass(),
                                    ref_core_ag.positions - ref_core_ag.center_of_mass())

    return np.arccos((R[0, 0] + R[1, 1] + R[2, 2] - 1) / 2) * 180 / np.pi


def get_normal_vector(leaflet):
    pos = leaflet.positions
    center = pos.mean(axis=0)
    posR = pos - center

    u, sigma, v = np.linalg.svd(posR)
    normal = v[2]
    normal = normal / np.linalg.norm(normal)
    if normal[2] < 0:
        normal = -normal
    return normal


def get_normal_vector_MIT(leafelet):
    return leafelet.principal_axes()[0]


def get_leaflets(u, cutoff=20):
    memb = u.select_atoms("resname POPE POPG and name P")

    valid = False

    while cutoff > 10:
        try:
            L = LeafletFinder(u, "resname POPE POPG and name P", cutoff=cutoff)
            group0 = L.groups(0)
            group1 = L.groups(1)
            assert len(memb) == len(group0 + group1)
            valid = True
            break
        except:
            print("Leaflet finder will be decreased from {} to {}".format(cutoff, cutoff-0.1))
            cutoff -= 0.1
            continue

    if not valid:
        raise ValueError("Could not find a suitable selection for lipids")

    if group0.center_of_mass()[2] > group1.center_of_mass()[2]:
        upper = group0
        lower = group1
    else:
        upper = group1
        lower = group0

    return upper, lower


def normal_to_spherical_projection(normal):

    x, y, z = normal

    theta = np.pi / 2 - np.arcsin(z)
    phi = np.arccos(x/np.sqrt(x**2 + y**2))
    phi = 2*np.pi - phi if y < 0 else phi
    return theta, phi


def collect(top, traj, filename=None, mutant=False):
    r"""Return a Dataset containing the timeseries of $\delta \phi$ and
    $\delta Z$ in both protomers. Runs at roughly 60 it/sec.

    Parameters:
    top (str): filename of topology definition
    traj (str): filename of trajectory file

    Keyword arguments:
    filename (str): output filename of data in netCDF format
                    Defaults to None.

    Returns:
    Dataset: xarray dataset of collective variables
    """
    if os.path.exists(str(filename)) and filename is not None:
        print("% already exists, please delete and run again" % filename)
        return -1

    from napalib.system.ref import OF_crystal_ref

    u = NapAUniverse(top, mutant=mutant)
    u.load_new(traj)
    nframes = len(u.trajectory)

    # allocate space
    dphi = np.zeros((nframes, 2))
    dz = np.zeros((nframes, 2))
    dtheta = np.zeros((nframes, 2))
    ddz = np.zeros((nframes, 2))
    dimer_z_memb = np.zeros((nframes, 2))
    core_z_memb = np.zeros((nframes, 2))
    membrane = np.zeros((nframes,))
    proj_phi = np.zeros((nframes,))
    proj_theta = np.zeros((nframes,))
    norm_data = np.zeros((nframes, 3))
    time = np.arange(0, nframes) * u.trajectory.dt

    membrane_ag = u.select_atoms("resname POPG POPE")

    c_traj_A = u.core_A.select_atoms("name CA")
    c_traj_B = u.core_B.select_atoms("name CA")
    d_traj_A = u.dimer_A.select_atoms("name CA")
    d_traj_B = u.dimer_B.select_atoms("name CA")

    c_ref = OF_crystal_ref.core_A.select_atoms("not bynum 4921 2221 and name CA")
    d_ref = OF_crystal_ref.dimer_A.select_atoms("not bynum 4921 2221 and name CA")

    upper, lower = get_leaflets(u)

    # Masrati CVs

    a1a = (c_traj_A + d_traj_A).select_atoms("resid 115 167 342 352")
    a2a = (c_traj_A + d_traj_A).select_atoms("resid 135 146 323 371")
    a3a = (c_traj_A + d_traj_A).select_atoms("resid 73  219 245 258")

    a1b = (c_traj_B + d_traj_B).select_atoms("resid 115 167 342 352")
    a2b = (c_traj_B + d_traj_B).select_atoms("resid 135 146 323 371")
    a3b = (c_traj_B + d_traj_B).select_atoms("resid 73  219 245 258")

    b1a = (c_traj_A + d_traj_A).select_atoms("resid 113-138 146-173 318-345")
    b2a = (c_traj_A + d_traj_A).select_atoms("resid 54-74 215-233 237-248 257-279")

    b1b = (c_traj_B + d_traj_B).select_atoms("resid 113-138 146-173 318-345")
    b2b = (c_traj_B + d_traj_B).select_atoms("resid 54-74 215-233 237-248 257-279")

    def com(x):
        return x.center_of_mass()

    for i, ts in tqdm(enumerate(u.trajectory), total=len(u.trajectory)):
        # need to collect z first because get_phi changes state...
        dz[i, 0] = get_z(c_traj_A, d_traj_A)
        dz[i, 1] = get_z(c_traj_B, d_traj_B)
        dimer_z_memb[i, 0] = get_z_memb(d_traj_A, membrane_ag)
        dimer_z_memb[i, 1] = get_z_memb(d_traj_B, membrane_ag)
        core_z_memb[i, 0] = get_z_memb(c_traj_A, membrane_ag)
        core_z_memb[i, 1] = get_z_memb(c_traj_B, membrane_ag)
        ddz[i, 0] = get_dz(com(b2a), com(b1a))
        ddz[i, 1] = get_dz(com(b2b), com(b1b))
        dtheta[i, 0] = get_theta(com(a1a), com(a2a), com(a3a))
        dtheta[i, 1] = get_theta(com(a1b), com(a2b), com(a3b))

        normal = get_normal_vector(upper)
        norm_data[i] = normal
        theta, phi = normal_to_spherical_projection(normal)

        proj_theta[i] = theta
        proj_phi[i] = phi

        membrane[i] = membrane_ag.center_of_mass()[2]

        dphi[i, 0] = get_phi(d_ref, c_ref, d_traj_A, c_traj_A)
        dphi[i, 1] = get_phi(d_ref, c_ref, d_traj_B, c_traj_B)

    data = xr.Dataset(
        {
            "dphi": (("time", "protomer"), dphi),
            "dtheta": (("time", "protomer"), dtheta),
            "ddz": (("time", "protomer"), ddz),
            "dz": (("time", "protomer"), dz),
            "dimer_z_memb": (("time", "protomer"), dimer_z_memb),
            "core_z_memb": (("time", "protomer"), core_z_memb),
            "memb_com": (("time"), membrane),
            "proj_phi": (("time"), proj_phi),
            "proj_theta": (("time"), proj_theta),
            "normal": (("time", "pos"), norm_data),
        },
        {
            "time": time,
            "protomer": ['A', 'B'],
            "pos": ['x', 'y', 'z'],
        },
    )

    if filename:
        data.to_netcdf(filename)

    return data
