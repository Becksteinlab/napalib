from napalib.system.universe import NapAUniverse

# load the topology
u = NapAUniverse("top.dms")

# load the trajectory, needs to be done separately for now
u.load_new("production.xtc")

# get atom groups from differing protomers and domains
core_A = u.core_A
core_B = u.core_B
dimer_A = u.dimer_A
dimer_B = u.dimer_B

# the often neglected helix 6, signified by neither domain
h6_A = u.neither_A
h6_B = u.neither_B

# THESE PROPERTIES ARE CALLING SELECTIONS, DO NOT USE IN A LOOP
# BIND THEM TO A VARIABLE AHEAD OF TIME

# Additional topology attributes include
# structure(s) -- "loop" or "helix"
# helix or helices -- names of helix "-1" "1" "2" ...
# domain -- "core" "dimer" "neither"
#
# These cannot be selected using select_atoms

# selection of helix 1 from both protomers
helix_1_both = u.residues[u.residues.helices == "-1"].atoms

# selection helix 1 from protomer A
helix_1_A = u.atoms[(u.atoms.helices == "-1") & (u.atoms.segids == "A")]
