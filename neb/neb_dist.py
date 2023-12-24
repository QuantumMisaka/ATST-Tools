# JamesMisaka in 2023-11-16
# read INIT and FINAL state and give the distance between each atoms
# part of ATST-Tools scripts

import os, sys
from ase.io import read, write
from ase.geometry import distance
from pathlib import Path

if __name__ == "__main__":
    msg = '''
Usage: 
    python neb_dist.py [init_image] [final_image]
Will return the distance of two image structure by Frobenius norm
'''
    if len(sys.argv) < 3:
        print(msg)
    else:
        init_infile = sys.argv[1]
        final_infile = sys.argv[2]
        init_Atoms = read(init_infile)
        final_Atoms = read(final_infile)
        # if constraints error, the constraint part need to be aligned
        # need more notice print-out
        distance = distance(init_Atoms, final_Atoms)
        print(f"==> Distance Between Two Image is {distance:.4f} Angstrom <==")