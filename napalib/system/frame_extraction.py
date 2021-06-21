from pathlib import Path
from napalib.system.traj import trajectories, Trajectory
from napalib.system.universe import NapAUniverse
from MDAnalysis import AtomGroup, Writer
from napalib.theme.IO import to_bendix_color_file, to_bendix_helix_file

def _write_coordinates(atoms: AtomGroup, filename: Path):
    """Write atoms to a file.

    Parameters
    ----------
    atoms: AtomGroup
        Target atoms for writing
    filename: Path
        Path for output file
    """
    with Writer(filename) as W:
        W.write(atoms)

def _extract_lipids(universe: NapAUniverse, outdir: Path):
    """Extract only POPE and POPG to a file.

    Parameters
    ----------
    universe: NapAUniverse
        Universe containing atomic positions.
    outdir: Path
        Base path
    """

    atoms = universe.select_atoms("resname POPE or resname POPG")
    _write_coordinates(atoms, outdir / "lipids.gro")

def _extract_protomer_A(universe: NapAUniverse, outdir: Path):
    """Extract only protomer A to a file.

    Parameters
    ----------
    universe: NapAUniverse
        Universe containing atomic positions.
    outdir: Path
        Base path
    """

    atoms = universe.atoms[universe.atoms.segids == "A"]
    _write_coordinates(atoms, outdir / "A.gro")

def _extract_protomer_B(universe: NapAUniverse, outdir: Path):
    """Extract only protomer B to a file.

    Parameters
    ----------
    universe: NapAUniverse
        Universe containing atomic positions.
    outdir: Path
        Base path
    """

    atoms = universe.atoms[universe.atoms.segids == "B"]
    _write_coordinates(atoms, outdir / "B.gro")

def _extract_solution(universe: NapAUniverse, outdir: Path):
    """Extract only solution to a file.

    Parameters
    ----------
    universe: NapAUniverse
        Universe containing atomic positions.
    outdir: Path
        Base path
    """

    atoms = universe.select_atoms("resname TIP3 or resname SOD or resname CLA")
    _write_coordinates(atoms, outdir / "solution.gro")

def extract_frame(trajectory: Trajectory, frame, outdir: Path):

    if not frame:
        return None

    outdir.mkdir(parents=True, exist_ok=True)

    u = NapAUniverse(trajectory.topology, mutant=True)
    u.load_new(trajectory.trajectory)

    u.trajectory[frame]

    _extract_lipids(u, outdir)
    _extract_protomer_A(u, outdir)
    _extract_protomer_B(u, outdir)
    _extract_solution(u, outdir)

    to_bendix_helix_file(outdir / "helices.txt")
    to_bendix_color_file(outdir / "colors.txt", dimmed=False)
    to_bendix_color_file(outdir / "colors_dimmed.txt", dimmed=True)

