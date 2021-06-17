from pathlib import Path
from napalib.system.universe import NapAUniverse


def to_bendix_helix_file(filename: Path):
    """Create a file containing the helix residue definitions that can be read using the Bendix plugin in VMD.
    The helix definitions are only defined for a single protomer.

    filename: Path
        Path of helix file
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

    filename: Path
        Path of color file
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
    """Return the selection schemes for mixed coloring of joining residues.

    :return: (str, str, str)
    """

    all_helices = NapAUniverse.core_helices + NapAUniverse.dimer_helices + NapAUniverse.cross_helix

    core = "resid "
    dimer = "resid "
    cross = "resid "

    def closest_domain(res):
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

    def in_helix(res):
        for lower, upper in all_helices:
            if res >= (lower+2) and res <= (upper - 2):
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