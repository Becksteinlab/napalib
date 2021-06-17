from pkg_resources import resource_filename
from napalib.system.universe import NapAUniverse
import os

IF_ref = resource_filename(__name__, 'if_ref.gro')
IF_crystal_ref = resource_filename(__name__, 'if_crystal_ref.pdb')
OF_crystal_ref = resource_filename(__name__, 'of_crystal_ref.pdb')
OCC_ref = resource_filename(__name__, 'occ_ref.gro')

IF_ref = NapAUniverse(IF_ref, mutant=True)
IF_crystal_ref = NapAUniverse(IF_crystal_ref, mutant=True)
OF_crystal_ref = NapAUniverse(OF_crystal_ref, mutant=True)
OCC_ref = NapAUniverse(OCC_ref, mutant=True)
