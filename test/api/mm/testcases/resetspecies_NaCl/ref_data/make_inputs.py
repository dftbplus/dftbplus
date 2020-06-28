"""  Import NaCl structure from DLPOLY and reorder
     Export the original and reordered structures in DFTB+ format
     using ASE.
     Alex Buccheri. 2020. alexanderbuccheri@googlemail.com
"""

import copy
import numpy as np
from ase.io import write, read

def sodium_chloride():

    nacl = read("NaCl512.config", format='dlp4')
    write("structure_1.gen", nacl)

    # Reorder atoms by reversing
    # (could also nudge positions if desired)
    nacl_reordered = copy.deepcopy(nacl)

    n_atoms = len(nacl.positions)
    index_map = np.arange(0, n_atoms)[::-1]

    # For API to work, need to ensure speciesName order is preserved i.e. Na Cl
    # This can be achieved by making sure the first atom in the list is Na
    index_map[0], index_map[1] = index_map[1], index_map[0]

    positions = []
    species = []
    for i in index_map:
        positions.append(nacl.positions[i])
        species.append(nacl.numbers[i])

    nacl_reordered.positions = positions
    nacl_reordered.numbers = species
    write("structure_2.gen", nacl_reordered)

    return

sodium_chloride() 
