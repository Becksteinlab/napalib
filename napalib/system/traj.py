from napalib.system.universe import NapAUniverse
import MDAnalysis as mda
from tqdm import tqdm
from pathlib import Path
from functools import reduce
from yaml import Loader, load
import os


class Trajchunk(object):

    def __init__(self, state, start, end, protomer):
        self.state = state
        self.start = start
        self.end = end
        self.protomer = protomer
        self.trajectory = None

    def to_dict(self):
        """Get dictionary representation of a Trajchunk.
        """
        export = {}
        export['state'] = self.state
        export['start'] = self.start
        export['end'] = self.end
        export['protomer'] = self.protomer
        return export

    @staticmethod
    def from_dict(data):
        """Create a Trajchunk from its dictionary representation.
        """
        state = str(data['state'])
        start = int(data['start'])
        end = int(data['end'])
        prot = str(data['protomer'])
        return Trajchunk(state, start, end, prot)


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

    def write_reduced(self, dirname, ions=True, water=False, lipids=False, stride=100, start=None, stop=None, verbose=True):
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

    def to_dict(self):
        """Get dictionary representation of a trajectory.
        """
        export = {}

        export['name'] = self._name
        export['top'] = self.topology
        export['traj'] = self.trajectory
        export['chunks'] = []

        for chunk in (self.chunks['A'] + self.chunks['B']):
            export['chunks'].append(chunk.to_dict())

        return export

    @classmethod
    def from_dict(cls, data):
        """Create Trajectory from its dictionary representation.
        """
        name = str(data['name'])
        topology = str(data['top'])
        trajectory = str(data['traj'])
        chunks = data['chunks']

        traj = cls(topology, trajectory, name)
        for chunk in chunks:
            traj.add_chunk(Trajchunk.from_dict(chunk))

        return traj

    def __str__(self):
        return self.name()

    def __repr__(self):
        return self.name()


class Anton2Trajectory(Trajectory):

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


def trajectories_from_file(filename, tag=None):
    """Load trajectory from YAML file.

    Parameters
    ----------
    filename : str
        Path to the YAML file containing trajectory definitions.

    Returns
    -------
    list
        List of trajectory objects
    """

    tagtable = {'general': Trajectory,
                'anton2': Anton2Trajectory,
                }

    with open(filename, 'r') as F:
        data = load(F, Loader)

    if data[0].get('tag', None):
        tag = data.pop(0)['tag']
    else:
        tag = 'general'

    return [tagtable[tag].from_dict(d) for d in data]


trajectories = list()

HOME_FILE = Path(os.environ['HOME']) / '.napalibtraj'
CWD_FILE = Path.cwd() / '.napalibtraj'

if HOME_FILE.exists():
    trajectories += trajectories_from_file(HOME_FILE)
elif CWD_FILE.exists():
    trajectories += trajectories_from_file(CWD_FILE)


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
