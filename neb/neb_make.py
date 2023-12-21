# JamesMisaka in 2023-11-02
# make NEB chain by ASE-NEB
# part of ATST-Tools scripts

import os, sys
# from ase.neb import NEB # old ase
from ase.mep.neb import NEB
from ase.io import read, write
from ase import Atoms
from ase.constraints import FixAtoms
from pathlib import Path


def set_fix_for_Atoms(atoms, fix_height=0, fix_dir=1,):
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

def get_neb_guess_chain(init_Atoms: Atoms, final_Atoms: Atoms, 
                  n_max: int, interpolate='idpp',
                  fix_height:float = None, fix_dir:int = 1, 
                  mag_ele:list = None, mag_num:list = None):
    """Get NEB chain by NEB method in ASE
    
    Args:
        init_Atoms (Atoms): initial state
        final_Atoms (Atoms): final state
        n_max (int): number of max images
        interpolate (string): interpolate chain path, 'linear' or 'idpp'
        fix_dir (int): fix atom in direction (0,1,2 for x,y,z), default 1
        fix_height (float): fix atom below height (fractional), default None
        mag_ele (list): elements for setting initial magmom for atoms, default None
        mag_num (list): number for setting initial magmom for atoms, default None
    """
    images = [init_Atoms]
    for i in range(n_max):
        image = init_Atoms.copy()
        # init and final image need to set mag and fix but not here
        if bool(fix_height):
            set_fix_for_Atoms(image, fix_height=fix_height, fix_dir=fix_dir)
        if bool(mag_ele):
            set_magmom_for_Atoms(image, mag_ele=mag_ele, mag_num=mag_num)
        images.append(image)
    images.append(final_Atoms)
    # use a simple NEB object to create traj file
    neb = NEB(images)
    # calculation and others should be set independently
    if interpolate in ['idpp','linear']:
        neb.interpolate(method=interpolate)
        # print-out guess information
    elif interpolate:
        print("---- Warning: interpolate method not supported, using default linear interpolate ----")
        neb.interpolate(method="linear") # using default
    return images


def nebmake(initial:Atoms=None, final:Atoms=None, ts_guess:Atoms=None,
                n_max:int=1, interpolate='idpp',  
                infile="input_guess_chain.traj", 
                outfile='init_neb_chain.traj', 
                fix_height:float = None, fix_dir:int = 1, 
                mag_ele:list = None, mag_num:list = None):
    """Make NEB Trajectory by ASE method, and print-out traj file, namely:
    1. images defining path from initial to final state
    2. neb object including the images and the implemented NEB method
    
    Args:
        initial (atoms): initial state object
        final (atoms): final state object
        n_max (int): number of max images
        interpolate (string): interpolate chain path, 'linear' or 'idpp'
        infile (string): input guess traj file, default name 'input_guess_chain.traj'
        outfile (string): output traj file, default name 'init_neb_chain.traj'
        fix_dir (int): fix atom in direction (0,1,2 for x,y,z), default 1
        fix_height (float): fix atom below height (fractional), default None
        mag_ele (list): elements for setting initial magmom for atoms, default None
        mag_num (list): number for setting initial magmom for atoms, default None
    """
    # make traj file from other input file is supported
    # for continuation
    if Path(infile).exists():
        print(f"----- Input Guess {infile} detected ! -----")
        print(f"----- Guess Trajectory just from {infile}  -----")
        if n_max > 0 and type(n_max) == int:
            images = read(f"{infile}@-{n_max + 2}:")
        else:
            raise ValueError("n_max should be a positive integer")
        # set constraint by atom height along some direction
        if bool(fix_height):
            print("---- Note: fix option is specified for input guess ----")
            for image in images[1:-1]:
                set_fix_for_Atoms(image, fix_height=fix_height, fix_dir=fix_dir)
        # set init-magmom for atom by element
        if bool(mag_ele):
            print("---- Note: init magmom option is specified for input guess ----")
            for image in images[1:-1]:
                set_magmom_for_Atoms(image, mag_ele=mag_ele, mag_num=mag_num)
        interpolate = None
        write(f'{outfile}', images, format='traj')
        # terminate nebmake for initial guess provided
        return
    elif (bool(initial) and bool(final)):
        if n_max > 0 and type(n_max) == int:
            # have TS_guess information
            if bool(ts_guess):
                print("----- Input TS guess detected and used ! -----")
                n_ts = n_max // 2 + 1
                n_init = n_ts - 1
                n_final = n_max - n_ts
                image_before_ts = get_neb_guess_chain(initial, ts_guess, n_init, interpolate,)
                image_after_ts = get_neb_guess_chain(ts_guess, final, n_final, interpolate,)
                images = image_before_ts + image_after_ts[1:]
            else:
            # not TS_guess information
                print("----- Get NEB guess chain only by initial and final state -----")
                images = get_neb_guess_chain(initial, final, n_max, interpolate,)
            print(f"--- Successfully make guessed image chain by {interpolate} method ! ---")
            write(f'{outfile}', images, format='traj')
            return images # for main function to read
        else:
            raise ValueError("n_max should be a positive integer")
    else:
        raise ValueError("initial and final state or input traj file should be provided")


if __name__ == "__main__":
    msg = '''
Usage: 
Default: 
    python neb_make.py [init_image] [final_image] [n_max] [optional]
[optional]:
    --fix [height]:[direction] : fix atom below height (fractional) in direction (0,1,2 for x,y,z), default None
    --mag [element1]:[magmom1],[element2]:[magmom2],... : set initial magmom for atoms of element, default None
    --format [format]: input image-files format, support abacus-out, vasp-out, extxyz, traj and others which have calculated property information, default None for automatically detect
    --ts [ts_guess]: use TS guess to make image chain, default None
Use existing guess and do continuation calculation
    python neb_make.py -i [init_guess_traj] [n_max] [optional]
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
        format = None
        ts_guess = None
        if "--fix" in sys.argv[4:]:
            fix_ind = sys.argv.index("--fix")
            height = float(sys.argv[fix_ind+1].split(":")[0])
            direction = int(sys.argv[fix_ind+1].split(":")[1])
        if "--mag" in sys.argv[4:]:
            mag_ind = sys.argv.index("--mag")
            mag_pair = sys.argv[mag_ind+1].split(",")
            mag_ele = [ele.split(":")[0] for ele in mag_pair]
            mag_num = [float(num.split(":")[1]) for num in mag_pair]
        if "--format" in sys.argv[4:]:
            format_ind = sys.argv.index("--format")
            format = sys.argv[format_ind+1]
        if "--ts" in sys.argv[4:]:
            ts_ind = sys.argv.index("--ts")
            try:
                ts_guess = read(sys.argv[ts_ind+1], format=format, )
            except:
                raise ValueError(f"Can't read TS guess file by specified format {format}")
        # set images
        if sys.argv[1] == "-i":
            infile : str = sys.argv[2]
            n_max = int(sys.argv[3])
            nebmake(infile=infile, n_max=n_max, interpolate=None, 
                    fix_height=height, fix_dir=direction, mag_ele=mag_ele, mag_num=mag_num)
        else:
            initial = sys.argv[1]
            final = sys.argv[2]
            init_Atoms = read(initial, format=format,)
            final_Atoms = read(final, format=format, )
            n_max = int(sys.argv[3])
            nebmake(init_Atoms, final_Atoms, ts_guess = ts_guess, 
                    n_max = n_max, interpolate='idpp', 
                    fix_height=height, fix_dir=direction, 
                    mag_ele=mag_ele, mag_num=mag_num)
    else:
        print(msg)



