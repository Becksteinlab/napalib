from napalib.system.universe import NapAUniverse
import MDAnalysis as mda
from tqdm import tqdm
from pathlib import Path
from functools import reduce


class Trajchunk(object):

    def __init__(self, state, start, end, protomer):
        self.state = state
        self.start = start
        self.end = end
        self.protomer = protomer
        self.trajectory = None


class Trajectory(object):

    def __init__(self, topology, trajectory, default_name):
        self._name = default_name
        self.topology = topology
        self.trajectory = trajectory
        self.chunks = {'A': [],
                       'B': []}

    def name(self, protomer=None):
        if protomer:
            return f"{self._name}_{protomer}"
        return f"{self._name}"

    def add_chunk(self, chunk):
        self.chunks[chunk.protomer.upper()].append(chunk)
        chunk.trajectory = self

    def universe(self):
        u = NapAUniverse(self.topology, mutant=True)
        u.load_new(self.trajectory)
        return u

    def write_reduced(self, dirname: Path, ions=True, water=False, lipids=False, stride=100, start=None, stop=None, verbose=True):
        """Write a reduced topology and trajectory into a directory.

        Parameters
        ----------
        dirname: Path
            Path to directory where files will be written.
        ions: bool
            Whether to include ions.
        water: bool
            Whether to include water.
        lipids: bool
            Whether to include lipids.
        stride: int
            Frame stride for write out.
        start: int
            Starting frame.
        stop: int
            Ending frame.
        verbose: bool
            Whether to use a progress bar.
        """
        u = self.universe()

        if not isinstance(dirname, Path):
            dirname = Path(dirname)

        dirname.mkdir(exist_ok=True, parents=True)
        top_file = dirname / "top.gro"
        traj_file = dirname / "traj.xtc"

        protein_sel = u.select_atoms("protein")
        lipids_sel = u.select_atoms("resname POPE POPG")
        water_sel = u.select_atoms("resname TIP3")
        ions_sel = u.select_atoms("resname SOD CLA")

        selections = [(protein_sel, True),
                      (lipids_sel, lipids),
                      (water_sel, water),
                      (ions_sel, ions)]

        included = list(filter(lambda x: x[1], selections))
        atoms = reduce(lambda x, y: x[0] + y[0], included)

        atoms.write(top_file)

        with mda.Writer(str(traj_file), atoms.n_atoms) as W:
            for ts in tqdm(u.trajectory[start:stop:stride], disable=(not verbose)):
                W.write(atoms)

    @property 
    def is_inward(self):
        _inward = 'inward' in self.name()
        _if = '_if_' in self.name()
        _ifm = '_ifm_' in self.name()
        return _inward or _if or _ifm

    @property
    def is_outward(self):
        _outward = 'outward' in self.name()
        _of = '_of_' in self.name()
        return _outward or _of

    @property
    def is_occluded(self):
        _OCC = '_OCC_' in self.name()
        return _OCC 

    @property
    def has_s2(self):
        return 's2' in self.name()

    @property
    def has_s1(self):
        return not self.has_s2

    @property
    def has_s4(self):
        return self.has_s1

    @property
    def is_310(self):
        return '310' in self.name()

    @property
    def is_358(self):
        return not self.is_310

    @property
    def if_s2_310(self):
        return self.is_inward and self.has_s2 and self.is_310


    @property
    def if_s4_310(self):
        return self.is_inward and self.has_s4 and self.is_310

    @property
    def of_s2_310(self):
        return self.is_outward and self.has_s2 and self.is_310


    @property
    def of_s4_310(self):
        return self.is_outward and self.has_s4 and self.is_310


    def __str__(self):
        return self.name()


    def __repr__(self):
        return self.name()


trajectories = list()

# anton_if_0
trajectories.append(Trajectory("/nfs/homes4/Projects/NapA/Anton2/traj/processed/direct/358/IF_WT/0/top.dms",
                               "/nfs/homes4/Projects/NapA/Anton2/traj/processed/direct/358/IF_WT/0/production.trr",
                               "a_if_0"))

trajectories[-1].add_chunk(Trajchunk("IF", 0, -1, 'A'))
trajectories[-1].add_chunk(Trajchunk("IF", 0, -1, 'B'))

# anton_occ_0
trajectories.append(Trajectory("/nfs/homes4/Projects/NapA/Anton2/traj/processed/direct/358/OCC_WT/0/top.dms",
                               "/nfs/homes4/Projects/NapA/Anton2/traj/processed/direct/358/OCC_WT/0/production.trr",
                               'a_occ_0'))

trajectories[-1].add_chunk(Trajchunk("OCC1", 0, 4860, 'A'))
trajectories[-1].add_chunk(Trajchunk("OCC1-IF", 4861, 4999, 'A'))
trajectories[-1].add_chunk(Trajchunk("IF", 5000, 17550, 'A'))
trajectories[-1].add_chunk(Trajchunk("IF-OCC2", 17551, 19010, 'A'))
trajectories[-1].add_chunk(Trajchunk("OCC2", 19011, -1, 'A'))

trajectories[-1].add_chunk(Trajchunk("OF", 0, -1, 'B'))

# anton_of_0
trajectories.append(Trajectory("/nfs/homes4/Projects/NapA/Anton2/traj/processed/direct/358/OF_WT/0/top.dms",
                               "/nfs/homes4/Projects/NapA/Anton2/traj/processed/direct/358/OF_WT/0/production.trr",
                               "a_of_0"))

trajectories[-1].add_chunk(Trajchunk("OF", 0, 14488, 'A'))
trajectories[-1].add_chunk(Trajchunk("OF-OCC", 14489, 15117, 'A'))
trajectories[-1].add_chunk(Trajchunk("OCC", 15118, 22192, 'A'))
trajectories[-1].add_chunk(Trajchunk("OCC-IF", 22193, 23241, 'A'))
trajectories[-1].add_chunk(Trajchunk("IF", 23242, -1, 'A'))

trajectories[-1].add_chunk(Trajchunk("OF", 0, -1, 'B'))

# anton_of_1
trajectories.append(Trajectory("/nfs/homes4/Projects/NapA/Anton2/traj/processed/direct/358/OF_WT/1/top.dms",
                               "/nfs/homes4/Projects/NapA/Anton2/traj/processed/direct/358/OF_WT/1/production.trr",
                               "a_of_1"))
trajectories[-1].add_chunk(Trajchunk("OF", 0, -1, 'A'))
trajectories[-1].add_chunk(Trajchunk("OF", 0, -1, 'B'))

# gromacs_if_1
trajectories.append(
    Trajectory("/nfs/homes4/Projects/NapA/Anton2/workspaces/ikenney/equilibrium_simulations/inward/01/md.tpr",
               "/nfs/homes4/Projects/NapA/Anton2/workspaces/ikenney/equilibrium_simulations/inward/01/production.xtc",
               "g_if_1"))

trajectories[-1].add_chunk(Trajchunk("IF", 0, -1, 'A'))
trajectories[-1].add_chunk(Trajchunk("IF", 0, -1, 'B'))

# gromacs_if_2
trajectories.append(
    Trajectory("/nfs/homes4/Projects/NapA/Anton2/workspaces/ikenney/equilibrium_simulations/inward/02/md.tpr",
               "/nfs/homes4/Projects/NapA/Anton2/workspaces/ikenney/equilibrium_simulations/inward/02/production.xtc",
               "g_if_2"))

trajectories[-1].add_chunk(Trajchunk("IF", 0, -1, 'A'))
trajectories[-1].add_chunk(Trajchunk("IF", 0, -1, 'B'))

# gromacs_if_3
trajectories.append(
    Trajectory("/nfs/homes4/Projects/NapA/Anton2/workspaces/ikenney/equilibrium_simulations/inward/03/md.tpr",
               "/nfs/homes4/Projects/NapA/Anton2/workspaces/ikenney/equilibrium_simulations/inward/03/production.xtc",
               "g_if_3"))

trajectories[-1].add_chunk(Trajchunk("IF", 0, -1, 'A'))
trajectories[-1].add_chunk(Trajchunk("IF", 0, 21500, 'B'))
trajectories[-1].add_chunk(Trajchunk("OF", 21501, -1, 'B'))

# gromacs_of_1
trajectories.append(
    Trajectory("/nfs/homes4/Projects/NapA/Anton2/workspaces/ikenney/equilibrium_simulations/outward/01/md.tpr",
               "/nfs/homes4/Projects/NapA/Anton2/workspaces/ikenney/equilibrium_simulations/outward/01/production.xtc",
               "g_of_1"))

trajectories[-1].add_chunk(Trajchunk("OF", 0, -1, 'A'))
trajectories[-1].add_chunk(Trajchunk("OF", 0, -1, 'B'))

# gromacs_of_2
trajectories.append(
    Trajectory("/nfs/homes4/Projects/NapA/Anton2/workspaces/ikenney/equilibrium_simulations/outward/02/md.tpr",
               "/nfs/homes4/Projects/NapA/Anton2/workspaces/ikenney/equilibrium_simulations/outward/02/production.xtc",
               "g_of_2"))

trajectories[-1].add_chunk(Trajchunk("OF", 0, -1, 'A'))
trajectories[-1].add_chunk(Trajchunk("OF", 0, -1, 'B'))

# gromacs_of_3
trajectories.append(
    Trajectory("/nfs/homes4/Projects/NapA/Anton2/workspaces/ikenney/equilibrium_simulations/outward/03/md.tpr",
               "/nfs/homes4/Projects/NapA/Anton2/workspaces/ikenney/equilibrium_simulations/outward/03/production.xtc",
               "g_of_3"))

trajectories[-1].add_chunk(Trajchunk("OF", 0, -1, 'A'))
trajectories[-1].add_chunk(Trajchunk("OF", 0, -1, 'B'))

# gromacs_occ_5

trajectories.append(Trajectory("/nfs/homes4/Projects/NapA/Anton2/traj/processed/direct/358/OF_WT/0/top.dms",
                               "/nfs/homes4/Projects/NapA/Anton2/workspaces/ikenney/chenou_traj/5/production.xtc",
                               "g_occ_5"))

trajectories[-1].add_chunk(Trajchunk("OCC", 0, 21500, 'A'))
trajectories[-1].add_chunk(Trajchunk("OCC-OF", 21501, 41999, 'A'))
trajectories[-1].add_chunk(Trajchunk("OF", 42000, -1, 'A'))

trajectories[-1].add_chunk(Trajchunk("OF", 0, -1, 'B'))

# gromacs_occ_*

for i in [1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20]:
    trajectories.append(Trajectory("/nfs/homes4/Projects/NapA/Anton2/traj/processed/direct/358/OF_WT/0/top.dms",
                                   f"/nfs/homes4/Projects/NapA/Anton2/workspaces/ikenney/chenou_traj/{i}/production.xtc",
                                   f"g_occ_{i}"))

    trajectories[-1].add_chunk(Trajchunk("OCC", 0, -1, 'A'))
    trajectories[-1].add_chunk(Trajchunk("OF", 0, -1, 'B'))


# low temp runs

def sort_parts(traj_list):
    print(traj_list)
    print(traj_list[0].split(".")[1][4:])
    traj_list.sort(key=lambda x: int(x.split(".")[1][4:]))
    return traj_list


for i in range(1, 11):
    for c in ["inward", "outward"]:
        conf = "if" if c == "inward" else "of"
        top = f"/nfs/homes4/Projects/NapA/Anton2/workspaces/ikenney/s2_310/{c}/production/{str(i).rjust(2, str(0))}/md.tpr"
        # _traj_list = glob(f"/nfs/homes4/Projects/NapA/Anton2/workspaces/ikenney/s2_310/inward/production/{str(i).rjust(2, str(0))}/*.xtc")
        traj = f"/nfs/homes4/Projects/NapA/Anton2/workspaces/ikenney/s2_310/{c}/production/{str(i).rjust(2, str(0))}/production.xtc"

        # if not _traj_list:
        #     continue

        # _traj_list = sort_parts(_traj_list)

        trajectories.append(Trajectory(top,
                                       traj,
                                       f"g_{c}_{i}_s2_310"))
        trajectories[-1].add_chunk(Trajchunk(conf.upper(), 0, -1, 'A'))
        trajectories[-1].add_chunk(Trajchunk(conf.upper(), 0, -1, 'B'))

for i in range(1, 4):
    for c in ['inward', 'outward']:
        conf = 'if' if c == 'inward' else 'of'
        top = f'/nfs/homes4/Projects/NapA/Anton2/workspaces/ikenney/s2_simulations/{c}/{str(i).rjust(2,str(0))}/gromacs/md.tpr'
        traj = f'/nfs/homes4/Projects/NapA/Anton2/workspaces/ikenney/s2_simulations/{c}/{str(i).rjust(2,str(0))}/gromacs/production.xtc'

        trajectories.append(Trajectory(top,
                                       traj,
                                       f'g_{c}_{i}_s2'))
        trajectories[-1].add_chunk(Trajchunk(conf.upper(), 0, -1, 'A'))
        trajectories[-1].add_chunk(Trajchunk(conf.upper(), 0, -1, 'B'))

# a_of_0_310
trajectories.append(Trajectory("/nfs/homes4/Projects/NapA/Anton2/traj/processed/direct/310/OF_WT/0/top.dms",
                               "/nfs/homes4/Projects/NapA/Anton2/traj/processed/direct/310/OF_WT/0/production.trr",
                               "a_of_0_310"))

trajectories[-1].add_chunk(Trajchunk("OF", 0, -1, 'A'))
trajectories[-1].add_chunk(Trajchunk("OF", 0, -1, 'B'))

# a_if_0_310
trajectories.append(Trajectory("/nfs/homes4/Projects/NapA/Anton2/traj/processed/direct/310/IF_WT/0/top.dms",
                               "/nfs/homes4/Projects/NapA/Anton2/traj/processed/direct/310/IF_WT/0/production.trr",
                               "a_if_0_310"))

trajectories[-1].add_chunk(Trajchunk("IF", 0, -1, 'A'))
trajectories[-1].add_chunk(Trajchunk("IF", 0, -1, 'B'))

# a_ifm_0_310
trajectories.append(Trajectory("/nfs/homes4/Projects/NapA/Anton2/traj/processed/direct/310/IF_MUT/0/top.dms",
                               "/nfs/homes4/Projects/NapA/Anton2/traj/processed/direct/310/IF_MUT/0/production.trr",
                               "a_ifm_0_310"))

trajectories[-1].add_chunk(Trajchunk("IF", 0, -1, 'A'))
trajectories[-1].add_chunk(Trajchunk("IF", 0, -1, 'B'))


def search_trajectories(traj_label):
    """Find a trajectory given the name of the trajectory.

    Parameters
    ----------
    traj_label: str
        Name of the trajectory.

    Returns
    -------
    Trajectory
        Trajectory object with the specified label
    """
    for trajectory in trajectories:
        if trajectory.name() == traj_label:
            return trajectory


def search_substates(substate):
    """Find all trajectory chunks that exist in the specified substate.

    Parameters
    ----------
    substate: str
        Substate to be matched.

    Returns
    -------
    [Trajchunk]
        All trajectory chunks in the specified substate.
    """
    states = []
    for trajectory in trajectories:
        for chunk in trajectory.chunks['A']:
            if chunk.state == substate:
                states.append(chunk)

        for chunk in trajectory.chunks['B']:
            if chunk.state == substate:
                states.append(chunk)

    return states


def is_310(traj):
    return '310' in traj.name()


def is_358(traj):
    return not is_310(traj)


def is_inward(traj):
    return 'inward' in traj.name() or 'if' in traj.name()


def is_outward(traj):
    return 'outward' in traj.name() or 'of' in traj.name()


def is_inward_310(traj):
    return is_310(traj) and is_inward(traj)


def is_outward_310(traj):
    return is_310(traj) and is_outward(traj)
