import MDAnalysis as mda
from napalib.system.universe import NapAUniverse
# from napalib.system.definitions import get_domain
import numpy as np

charge_groups = {
    "LYS": ["NZ"],
    "GLU": ["OE1", "OE2"],
    "ASP": ["OD1", "OD2"],
    "ARG": ["NH2"]
}

positive = ["LYS", "ARG"]
negative = ["GLU", "ASP"]


class Residue(object):
    def __init__(self, resname, resid):
        if not resname in charge_groups.keys():
            print("Could not find %s in valid residues" % resname)
            raise ValueError
        self.resname = resname
        self.resid = resid
        self._atoms = []
        self.positive = True if self.resname in positive else False
        self.negative = not self.positive

    @staticmethod
    def from_label(label):
        resname = label[0:3]
        resid = int(label[3:])
        return Residue(resname, resid)

    @property
    def atoms(self):
        return self._atoms

    @atoms.setter
    def atoms(self, atoms):
        if len(atoms) != len(charge_groups[self.resname]):
            raise ValueError
        self._atoms = atoms

    def get_domain(self):
        return NapAUniverse._get_domain(self.resid)

    def __eq__(self, other):
        return self.resname == other.resname and self.resid == other.resid

    def __repr__(self):
        return self.resname + str(self.resid)

    def __str__(self):
        return self.__repr__()


class Pair(object):
    """
    
    """

    def __init__(self, R1, R2):
        assert type(R1) is Residue
        assert type(R2) is Residue

        self.positive = R1 if R1.resname in positive else R2
        self.negative = R2 if self.positive == R1 else R1

        if self.positive.resname == self.negative.resname:
            raise ValueError("You cannot have a salt bridge between two identical residues")

    def calculate_separation(self):
        return min([np.linalg.norm(p.position - n.position) for p in self.positive.atoms for n in self.negative.atoms])

    def interdomain(self):
        return not (self.positive.get_domain() == self.negative.get_domain())

    def get_time_series(self, da):
        return da.sel(pos=str(self.positive), neg=str(self.negative))

    def __repr__(self):
        return self.positive.__repr__() + " <-> " + self.negative.__repr__()

    def __eq__(self, RHS):
        return (RHS.positive == self.positive) and (RHS.negative == self.negative)

    @property
    def mapping(self):
        return {'pos': str(self.positive), 'neg': str(self.negative)}


class PairList():
    def __init__(self, pairs):
        self.pairs = pairs

    @property
    def positive_labels(self):
        return list(set([str(p.positive) for p in self.pairs]))

    @property
    def negative_labels(self):
        return list(set([str(p.negative) for p in self.pairs]))

    @property
    def mapping(self):
        pos = self.positive_labels
        neg = self.negative_labels
        return {'pos': pos, 'neg': neg}

    def existence_matrix(self, da, cutoff=5):
        d = da.sel(**self.mapping)
        return (d < cutoff).sum(axis=0) / da.shape[0]


def get_charged_residues(universe):
    """Returns a dictionary with the residue IDs from each protomer given the charge
    
    """
    results = {"A": [], "B": []}

    ag_A = universe.select_atoms("segid A")
    ag_B = universe.select_atoms("segid B")

    for residue in ag_A.residues + ag_B.residues:
        residue_charge = sum(residue.atoms.charges)
        if abs(residue_charge) > 0.2 and residue.resname in charge_groups.keys():
            R = Residue(residue.resname, residue.resid)
            R.atoms = sum([residue.atoms.select_atoms("name " + name) for name in charge_groups[R.resname]])
            results[residue.segid].append(R)
    return results


def collection_scheme(residue_dict):
    results = {'A': [], 'B': []}

    group_A = residue_dict['A']
    group_B = residue_dict['B']

    for p in filter(lambda x: x.positive, group_A):
        for n in filter(lambda x: x.negative, group_A):
            results['A'].append(Pair(p, n))

    for p in filter(lambda x: x.positive, group_B):
        for n in filter(lambda x: x.negative, group_B):
            results['B'].append(Pair(p, n))

    return results
