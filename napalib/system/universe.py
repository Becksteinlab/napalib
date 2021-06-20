import MDAnalysis as mda
import numpy as np

from MDAnalysis.core import topologyattrs


class SecondaryStructure(topologyattrs.ResidueAttr):
    attrname = "structures"
    singular = "structure"
    dtype = object

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.array(['' for _ in range(nr)], dtype=object)


class Helix(topologyattrs.ResidueAttr):
    attrname = "helices"
    singular = "helix"
    dtype = object

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.array(['' for _ in range(nr)], dtype=object)


class Domain(topologyattrs.ResidueAttr):
    attrname = "domains"
    singular = "domain"
    dtype = object

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.array(['' for _ in range(nr)], dtype=object)


class NapAUniverse(mda.Universe):
    dimer_helices = [(3, 27),
                     (33, 43),
                     (54, 74),
                     (214, 233),
                     (236, 247),
                     (257, 279)]

    dimer_helix_names = ["-1",
                         "1",
                         "2",
                         "7",
                         "8",
                         "9"]

    core_helices = [(86, 109),
                    (112, 125),
                    (127, 138),
                    (143, 173),
                    (287, 312),
                    (317, 328),
                    (331, 346),
                    (350, 385)]

    core_helix_names = ["3",
                        "4a",
                        "4b",
                        "5",
                        "10",
                        "11a",
                        "11b",
                        "12"]

    cross_helix = [(177, 201)]
    cross_helix_name = ["6"]

    @staticmethod
    def _selection_string_from_ranges(ranges):
        ids = [str(i) + "-" + str(j) for i, j in ranges]
        return "resid " + ' '.join(ids)

    @staticmethod
    def core_selection_string():
        helices = NapAUniverse.core_helices
        return NapAUniverse._selection_string_from_ranges(helices)

    @staticmethod
    def dimer_selection_string():
        helices = NapAUniverse.dimer_helices
        return NapAUniverse._selection_string_from_ranges(helices)

    @staticmethod
    def _get_domain(resid):
        for lower, upper in NapAUniverse.core_helices:
            if resid >= lower or resid <= upper:
                return "core"

        for lower, upper in NapAUniverse.dimer_helices:
            if resid >= lower or resid <= upper:
                return "dimer"

        return "neither"

    @staticmethod
    def _expand(i, j):
        return [j for j in range(i, j + 1)]

    def __init__(self, topology_file, mutant=False):
        super().__init__(topology_file)
        self.format = topology_file.split('.')[-1]
        self.mutant = mutant
        self._fix_topology()
        self._bind_secondary_structure()
        self._bind_helix_names()
        self._bind_domain()

    def _fix_topology(self):
        """Make the topology consistent between multiple filetypes.
        """

        protein = self.select_atoms("protein")
        nres = len(protein.residues)

        if self.format == "tpr":
            # This is typically the condition that we see with TPRs.
            # This code will needed to be edited if this ever fails
            assert len(protein.residues) - 1 == protein.residues[-1].resnum

            # TPRs have extended resids instead of the intended wrapping 
            # effect for the two protomers. 
            nres = len(protein.residues)
            assert protein.resids[0] == 0
            assert protein.resids[-1] == nres - 1
            # split in two even parts and offset by 3
            ids = list(protein.residues.resids % int(nres / 2) + 3)
            protein.residues.resids = ids

            protein.residues[:ids[-1] - 1].segments.segids = "A"
            protein.residues[ids[-1] - 1:].segments.segids = "B"

        elif self.format == "dms":
            # The most common problem with dms files is that the 
            # segids have not been changed to reflect the chains
            # also check that the development branch is installed

            assert list(np.unique(protein.residues.segids)) == ["A", "B"], \
                "Are you using the development version of MDAnalysis?"

        elif self.format == "gro":
            # GRO files need to have the segids added explicitly
            a = self.add_Segment(segid='A')
            b = self.add_Segment(segid='B')
            protein.residues[:int(nres / 2)].segments = a
            protein.residues[int(nres / 2):].segments = b

        elif self.format == "pdb":
            # PDBs are only really going to be the crystal structures...
            pass

        if self.format != "pdb":
            a = self.select_atoms("segid A")
            b = self.select_atoms("segid B")

            for i, j in zip(a.residues, b.residues):
                if not self.mutant:
                    assert i.resname == j.resname
                assert i.resid == j.resid

            assert len(a.residues) == len(b.residues)
        else:
            pass

    def _bind_secondary_structure(self):
        """Assign either loop or helix to each residue.

        """
        values = ["" for _ in range(len(self.residues))]
        self.add_TopologyAttr(SecondaryStructure(values))

        prot = self.select_atoms("protein")
        prot.residues.structures = "loop"
        helices = self.dimer_helices + self.core_helices + self.cross_helix

        targets = np.concatenate([NapAUniverse._expand(i, j) for i, j in helices])
        for t in targets:
            prot.residues[prot.residues.resids == t].structures = "helix"

    def _bind_helix_names(self):
        """Assign helix names to each residue.

        """
        values = ["" for _ in range(len(self.residues))]
        self.add_TopologyAttr(Helix(values))

        prot = self.select_atoms("protein")

        helices = [zip(self.dimer_helices, self.dimer_helix_names),
                   zip(self.core_helices, self.core_helix_names),
                   zip(self.cross_helix, self.cross_helix_name)]

        for hdefs in helices:
            for hrange, name in hdefs:
                targets = NapAUniverse._expand(hrange[0], hrange[1])
                targets = np.array(targets, dtype=object)
                for t in targets:
                    prot.residues[prot.residues.resids == t].helices = name

    def _bind_domain(self):
        """Assign the domain to each residue.

        """
        values = ["" for _ in range(len(self.residues))]
        self.add_TopologyAttr(Domain(values))

        protomer = self.select_atoms("protein")

        for name in self.dimer_helix_names:
            protomer.residues[protomer.residues.helices == name].domains = "dimer"

        for name in self.core_helix_names:
            protomer.residues[protomer.residues.helices == name].domains = "core"

        for name in self.cross_helix_name:
            protomer.residues[protomer.residues.helices == name].domains = "neither"

    def domain_select(self, protomer, domain):
        """Select an AtomGroup from a given domain in a given protomer.

        Parameters
        ----------
        protomer: str
            Protomer to select the atom group from. This can be either "A" or "B"
        domain: str
            Domain name. This can be either "core", "dimer", or "neither".

        Returns
        -------
        AtomGroup
            Atom group from selection
        """
        if self.format == "pdb":
            protomer = "A"
        ag = self.select_atoms("protein and segid " + protomer)
        ag = ag.atoms[ag.atoms.domains == domain]
        return ag

    @property
    def core_A(self):
        return self.domain_select("A", "core")

    @property
    def core_B(self):
        return self.domain_select("B", "core")

    @property
    def dimer_A(self):
        return self.domain_select("A", "dimer")

    @property
    def dimer_B(self):
        return self.domain_select("B", "dimer")

    @property
    def neither_A(self):
        return self.domain_select("A", "neither")

    @property
    def neither_B(self):
        return self.domain_select("B", "neither")
