# JamesMisaka in 2023-11-02
# make NEB chain by ASE-NEB
# part of ASE-NEB-ABACUS scripts

import os, sys
# from ase.neb import NEB # old ase
from ase.mep.neb import NEB
# from ase.autoneb import AutoNEB
from ase.io import read, write
from ase.constraints import FixAtoms
from pathlib import Path


def set_Atoms(infile, format=None, fix_height=0.5, fix_dir=1, 
              mag_ele=[], mag_num=[]):
    """Set Atoms Object by read from file
    
    Args:
        infile: (string) read-in Atoms file
        format (None or string) infile format, None for auto-detected
    """
    atoms = read(infile, format=format)
    # maybe we should set constaints independently in here
    if not fix_height:
        mask = atoms.get_scaled_positions()[:,fix_dir] < fix_height
        fix = FixAtoms(mask=mask)
        atoms.set_constraint(fix)
    # maybe we should set init magmom independently in here
    if not mag_ele:
        init_magmom = atoms.get_initial_magnetic_moments()
        if len(mag_ele) != len(mag_num):
            raise SyntaxWarning("mag_ele and mag_num have different length")
        for mag_pair in zip(mag_ele, mag_num):
            ele_ind = [atom.index for atom in atoms if atom.symbol == mag_pair[0]]
            init_magmom[ele_ind] = mag_num
        atoms.set_initial_magnetic_moments(init_magmom)
    return atoms


def nebmake(initial='', final='', n_max=8, interpolate='idpp',  
                outfile='init_neb_chain.traj', ):
    """Make NEB Trajectory by ASE method, and print-out traj file, namely:
    1. images defining path from initial to final state
    2. neb object including the images and the implemented NEB method
    
    Args:
        initial (atoms): initial state object
        final (atoms): final state object
        n_max (int): number of max images
        interpolate (string): interpolate chain path, 'linear' or 'idpp' or None
    """
    # outfile should not exists, however contination is supported
    if Path(outfile).exists():
        print(f"----- Target traj file {outfile} exist ! -----")
        print(f"----- Guess Trajectory just from {outfile}  -----")
        images = read(f"{outfile}@-{n_max + 2}:")
        interpolate = None
        return neb 
        # terminate nebmake for initial guess provided
    else:
        images = [initial]
        for i in range(n_max):
            image = initial.copy()
            images.append(image)
        images.append(final)
        # use a simple NEB object to create traj file
        neb = NEB(images)
        # calculation and others should be set independently
        if interpolate in ['idpp','linear']:
            neb.interpolate(method=interpolate)
            # print-out guess information
        elif interpolate:
            print("---- Warning: interpolate method not supported, using default linear interpolate ----")
            neb.interpolate(method="linear") # using default
        print(f"--- Successfully make guessed image chain by {interpolate} method ! ---")
        write(f'{outfile}', images, format='traj')
        return neb # for main function to read


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python neb_make.py [init_image] [final_image] [n_max] ")
    elif len(sys.argv) >= 4:
        initial = sys.argv[1]
        final = sys.argv[2]
        nmax = int(sys.argv[3])
        init_Atoms = set_Atoms(initial, format="abacus-out")
        final_Atoms = set_Atoms(final, format="abacus-out")
        nebmake(init_Atoms, final_Atoms, nmax)




