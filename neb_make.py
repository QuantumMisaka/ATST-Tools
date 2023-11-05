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


def set_Atoms(infile, format=None, fix_height=0.0, fix_dir=1, 
              mag_ele=[], mag_num=[]):
    """Set Atoms Object by read from file
    
    Args:
        infile: (string) read-in Atoms file
        format (None or string) infile format, None for auto-detected
    """
    atoms = read(infile, format=format)
    # maybe we should set constaints independently in here
    if fix_height:
        if fix_dir not in [0, 1, 2]:
            raise ValueError("fix_dir should be 0, 1 or 2 for x, y or z")
        mask = atoms.get_scaled_positions()[:, fix_dir] < fix_height
        fix = FixAtoms(mask=mask)
        atoms.set_constraint(fix)
    # maybe we should set init magmom independently in here
    if mag_ele:
        init_magmom = atoms.get_initial_magnetic_moments()
        if len(mag_ele) != len(mag_num):
            raise SyntaxWarning("mag_ele and mag_num have different length")
        for mag_pair in zip(mag_ele, mag_num):
            ele_ind = [atom.index for atom in atoms if atom.symbol == mag_pair[0]]
            init_magmom[ele_ind] = mag_num
        atoms.set_initial_magnetic_moments(init_magmom)
    return atoms


def nebmake(initial='', final='', n_max=8, interpolate='idpp',  
                infile="input_guess_chain.traj", outfile='init_neb_chain.traj', ):
    """Make NEB Trajectory by ASE method, and print-out traj file, namely:
    1. images defining path from initial to final state
    2. neb object including the images and the implemented NEB method
    
    Args:
        initial (atoms): initial state object
        final (atoms): final state object
        n_max (int): number of max images
        interpolate (string): interpolate chain path, 'linear' or 'idpp' or None
        infile (string): input guess traj file, default name 'input_guess_chain.traj'
        outfile (string): output traj file, default name 'init_neb_chain.traj'
    """
    # make traj file from other input file is supported
    if Path(infile).exists():
        print(f"----- Input Guess {infile} detected ! -----")
        print(f"----- Guess Trajectory just from {infile}  -----")
        images = read(f"{infile}@-{n_max + 2}:")
        interpolate = None
        return neb 
        # terminate nebmake for initial guess provided
    else:
        images = [initial]
        if nmax == 0:
            # always for autoneb, just do format transfer
            print("---- Warning: n_max = 0, traj file only contain initial and final images. ---- \n ---- If you are using AutoNEB method, just ignore it ----")
            images.append(final)
            write(f'{outfile}', images, format='traj')
        elif nmax > 0 and type(nmax) == int:
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
        else:
            raise ValueError("n_max should be a non-negative integer")


if __name__ == "__main__":
    msg = '''
Usage: 
Default: 
    python neb_make.py [init_image] [final_image] [n_max] [optional]
[optional]:
    --fix [height]:[direction] : fix atom below height (fractional) in direction (0,1,2 for x,y,z)
    --mag [element1]:[magmom1],[element2]:[magmom2],... : set initial magmom for atoms of element
Use existing guess: 
    python neb_make.py -i [init_guess_traj]
'''
    if len(sys.argv) < 3:
        print(msg)
        exit()
    elif len(sys.argv) == 3:
        if sys.argv[1] == "-i":
            infile = sys.argv[2]
            nebmake(infile=infile)
        else:
            print(msg)
            exit()
    elif len(sys.argv) >= 4:
        initial = sys.argv[1]
        final = sys.argv[2]
        nmax = int(sys.argv[3])
        height = 0.0
        direction = 1
        mag_ele = []
        mag_num = []
        if "--fix" in sys.argv[4:]:
            fix_ind = sys.argv.index("--fix")
            height = float(sys.argv[fix_ind+1].split(":")[0])
            direction = int(sys.argv[fix_ind+1].split(":")[1])
        if "--mag" in sys.argv[4:]:
            mag_ind = sys.argv.index("--mag")
            mag_pair = sys.argv[mag_ind+1].split(",")
            mag_ele = [ele.split(":")[0] for ele in mag_pair]
            mag_num = [float(num.split(":")[1]) for num in mag_pair]
        init_Atoms = set_Atoms(initial, format="abacus-out", 
                            fix_height=height, fix_dir=direction, mag_ele=mag_ele, mag_num=mag_num)
        final_Atoms = set_Atoms(final, format="abacus-out", 
                            fix_height=height, fix_dir=direction, mag_ele=mag_ele, mag_num=mag_num)
        nebmake(init_Atoms, final_Atoms, nmax)
    else:
        print(msg)



