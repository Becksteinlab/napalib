from pathlib import Path
from napalib.system.universe import NapAUniverse


def to_bendix_helix_file(filename: Path):
    """Create a file containing the helix residue definitions that can be read using the Bendix plugin in VMD.
    The helix definitions are only defined for a single protomer.

    Parameters
    ----------
    filename: Path
        Path of output helix file.
    """

    helices = NapAUniverse.core_helices + NapAUniverse.dimer_helices + NapAUniverse.cross_helix
    helices.sort(key=lambda x: x[0])

    contents = ""

    for start, stop in helices:
        contents += f"{start} {stop} "

    with open(filename, 'w') as F:
        F.write(contents)


def to_bendix_color_file(filename: Path, dimmed=False):
    """Create a file containing helix color definitions that can be read using the Bendix plugin in VMD.

    The colors need to be overwritten in VMD if you're not happy with a violet, magenta, and orange color-scheme.

    Parameters
    ----------
    filename: Path
        Path of output color file.
    dimmed: bool
        Whether or not to use the dimmed variant of the color scheme.
    """

    _colors = {'core': "violet2",
               'dimer': "magenta2",
               'cross': "orange2"}

    _colors_dimmed = {'core': "violet",
                      'dimer': "magenta",
                      'cross': "orange"}

    colors = _colors_dimmed if dimmed else _colors

    core_helices = [(i, colors['core']) for i in NapAUniverse.core_helices]
    dimer_helices = [(i, colors['dimer']) for i in NapAUniverse.dimer_helices]
    cross_helix = [(i, colors['cross']) for i in NapAUniverse.cross_helix]

    helices = core_helices + dimer_helices + cross_helix
    helices.sort(key=lambda x: x[0][0])  # sort by the first residue in a helix

    # extract the colors for the sorted list and join
    contents = " ".join(map(lambda x: x[1], helices))

    with open(filename, 'w') as F:
        F.write(contents)

def joining_residue_selections():
    """Assign non-helix secondary structure residues to a domain depending on the residues distance to those domains.

    Returns
    -------
    (str, str, str)
        VMD readable selection strings for the core, dimer, and cross-domain loops.
    """

    all_helices = NapAUniverse.core_helices + NapAUniverse.dimer_helices + NapAUniverse.cross_helix

    core = "resid "
    dimer = "resid "
    cross = "resid "

    def closest_domain(res):
        """Determine the closest domain to a residue.

        Parameters
        ----------
        res: int
            Resid

        Returns
        -------
        str
            Domain closest to the residue with resid res.

        """
        closest_distance = 385
        _closest_domain = None

        for helix in NapAUniverse.core_helices:
            for side in helix:
                distance = abs(side - res)
                if distance < closest_distance:
                    closest_distance = distance
                    _closest_domain = "core"

        for helix in NapAUniverse.dimer_helices:
            for side in helix:
                distance = abs(side - res)
                if distance < closest_distance:
                    closest_distance = distance
                    _closest_domain = "dimer"

        for helix in NapAUniverse.cross_helix:
            for side in helix:
                distance = abs(side - res)
                if distance < closest_distance:
                    closest_distance = distance
                    _closest_domain = "cross"

        return _closest_domain

    def in_helix(res, depth=2):
        for lower, upper in all_helices:
            if res >= (lower+depth) and res <= (upper - depth):
                return True
        return False

    for residue in filter(lambda x: not in_helix(x), range(3, 385)):
        domain = closest_domain(residue)

        if domain == "core":
            core += f"{residue} "

        elif domain == "dimer":
            dimer += f"{residue} "

        elif domain == "cross":
            cross += f"{residue} "

    return core, dimer, cross