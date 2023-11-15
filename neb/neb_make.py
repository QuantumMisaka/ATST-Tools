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


def set_fix_for_Atoms(atoms, fix_height=0.0, fix_dir=1,):
    # maybe we should set constaints independently in here
    if fix_dir not in [0, 1, 2]:
        raise ValueError("fix_dir should be 0, 1 or 2 for x, y or z")
    direction = {0: "x", 1: "y", 2: "z"}
    mask = atoms.get_scaled_positions()[:, fix_dir] <= fix_height
    fix = FixAtoms(mask=mask)
    print(f"---- Fix Atoms below {fix_height} in direction {direction[fix_dir]} ----")
    atoms.set_constraint(fix)



def set_magmom_for_Atoms(atoms, mag_ele=[], mag_num=[]):
    """Set Atoms Object magmom by element
    
    Args:
        atoms: (atoms) Atoms object
        mag_ele (list): element list
        mag_num (list): magmom list
    """
    # init magmom can only be set to intermediate images
    init_magmom = atoms.get_initial_magnetic_moments()
    if len(mag_ele) != len(mag_num):
        raise SyntaxWarning("mag_ele and mag_num have different length")
    for mag_pair in zip(mag_ele, mag_num):
        ele_ind = [atom.index for atom in atoms if atom.symbol == mag_pair[0]]
        init_magmom[ele_ind] = mag_pair[1]
    atoms.set_initial_magnetic_moments(init_magmom)
    print(f"---- Set initial magmom for {mag_ele} to {mag_num} ----")


def nebmake(initial, final, n_max, interpolate='idpp',  
                infile="input_guess_chain.traj", outfile='init_neb_chain.traj', 
                fix_height=0.0, fix_dir=1, mag_ele=[], mag_num=[]):
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
    # for continuation
    if Path(infile).exists():
        print(f"----- Input Guess {infile} detected ! -----")
        print(f"----- Guess Trajectory just from {infile}  -----")
        images = read(f"{infile}@-{n_max + 2}:")
        # set constraint by atom height along some direction
        if fix_height:
            for image in images:
                set_fix_for_Atoms(image, fix_height=fix_height, fix_dir=fix_dir)
        else:
            print("---- Warning: no fix height provided, no constraint set in ASE ----")
        # set init-magmom for atom by element
        if mag_ele:
            for image in images:
                set_magmom_for_Atoms(image, mag_ele=mag_ele, mag_num=mag_num)
        else:
            print("---- Warning: no element for initial magmom height provided, no magmom set in ASE ----")
        interpolate = None
        write(f'{outfile}', images, format='traj')
        # terminate nebmake for initial guess provided
        return
    else:
        images = [initial]
        if n_max == 0:
            # always for autoneb, just do format transfer
            print("---- Warning: n_max = 0, traj file only contain initial and final images.")
            images.append(final)
            write(f'{outfile}', images, format='traj')
        elif n_max > 0 and type(n_max) == int:
            for i in range(n_max):
                image = initial.copy()
                set_fix_for_Atoms(image, fix_height=fix_height, fix_dir=fix_dir)
                set_magmom_for_Atoms(image, mag_ele=mag_ele, mag_num=mag_num)
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
    python neb_make.py -i [input_guess_traj] [n_max]
Notice: n_max is related to the number of interpolated image, not include initial and final image
'''
    if len(sys.argv) < 4:
        print(msg)
        exit()
    elif len(sys.argv) >= 4:
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
        # set images
        if sys.argv[1] == "-i":
            infile = sys.argv[2]
            n_max = int(sys.argv[3])
            nebmake(infile=infile, n_max=n_max, interpolate=None, 
                    fix_height=height, fix_dir=direction, mag_ele=mag_ele, mag_num=mag_num)
        else:
            initial = sys.argv[1]
            final = sys.argv[2]
            init_Atoms = read(initial, format="abacus-out",)
            final_Atoms = read(final, format="abacus-out", )
            n_max = int(sys.argv[3])
            nebmake(init_Atoms, final_Atoms, n_max, interpolate='idpp', 
                    fix_height=height, fix_dir=direction, mag_ele=mag_ele, mag_num=mag_num)
    else:
        print(msg)



