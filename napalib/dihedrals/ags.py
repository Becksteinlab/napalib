from MDAnalysis.core.groups import AtomGroup


def phi_pairs(RG):
    return list(zip(RG[1:-1], RG[:-2]))


def psi_pairs(RG):
    return list(zip(RG[1:-1], RG[2:]))


def get_chi_ag(residue):
    if residue.resname in ("ALA", "GLY"):
        return None
    backbone = residue.atoms.select_atoms("backbone")
    n = list(backbone.select_atoms("name N"))
    ca = list(backbone.select_atoms("name CA"))
    cb = list(residue.atoms.select_atoms("name CB"))
    cg_str = "CG"
    if residue.resname in ("VAL", "ILE"):
        cg_str = "CG1"
    if residue.resname == "THR":
        cg_str = "OG1"
    if residue.resname == "SER":
        cg_str = "OG"
    if residue.resname == "CYS":
        cg_str = "SG"
    cg = list(residue.atoms.select_atoms("name {0}".format(cg_str)))
    ag = AtomGroup(n + ca + cb + cg)
    if len(ag) != 4:
        print(residue)
        raise ValueError("wrong number in selection")
    return ag


def get_phi_ag(pair):
    previous = pair[1].atoms
    current = pair[0].atoms
    prev_c = list(previous.select_atoms("backbone and name C"))
    others = list(current.select_atoms("backbone and not name O"))
    return AtomGroup(prev_c + others)


def get_psi_ag(pair):
    following = pair[1].atoms
    current = pair[0].atoms
    current_atoms = list(current.select_atoms("backbone and not name O"))
    following_N = list(following.select_atoms("backbone and name N"))
    return AtomGroup(current_atoms + following_N)
