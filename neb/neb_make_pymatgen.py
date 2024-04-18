import argparse
from ase.io import read, write
from ase import Atoms
from ase.constraints import FixAtoms
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.diffusion.neb.pathfinder import IDPPSolver
from ase.calculators.singlepoint import SinglePointCalculator

'''
Author: mosey
Last update: 2024-03-28

Dependencies:
1. Install pymatgen: pip install pymatgen (or conda install -c conda-forge pymatgen)
2. Install pymatgen-analysis-diffusion: pip install pymatgen-analysis-diffusion
3. Install ase-abacus interfaces
'''

def set_fix_for_Atoms(atoms: Atoms, fix_height: float=0, fix_dir: int=1,):
    # maybe we should set constaints independently in here
    if fix_dir not in [0, 1, 2]:
        raise ValueError("fix_dir should be 0, 1 or 2 for x, y or z")
    direction = {0: "x", 1: "y", 2: "z"}
    mask = atoms.get_scaled_positions()[:, fix_dir] <= fix_height
    fix = FixAtoms(mask=mask)
    print(f"Fix Atoms below {fix_height} in direction {direction[fix_dir]}")
    atoms.set_constraint(fix)

def set_magmom_for_Atoms(atoms: Atoms, mag_ele: list=[], mag_num: list=[]):
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
    print(f"Set initial magmom for {mag_ele} to {mag_num}")

def generate(method:str, n_images:int, is_file:str, fs_file:str, 
             output_file:str,sort_tol:float,format:str,
             fix_height: float=None, fix_dir: int=None,
             mag_ele: list=None, mag_num: list=None):

    # read files
    print(f'Reading files: {is_file} and {fs_file}')
    is_atom = read(is_file,format=format)
    fs_atom = read(fs_file,format=format)
    is_e = is_atom.get_potential_energy()
    fs_e = fs_atom.get_potential_energy()
    is_f = is_atom.get_forces()
    fs_f = fs_atom.get_forces()
    is_pmg = AseAtomsAdaptor.get_structure(is_atom)
    fs_pmg = AseAtomsAdaptor.get_structure(fs_atom)

    # generate path
    print(f'Generating path, number of images: {n_images}, sort_tol: {sort_tol}')

    # optimize path
    print(f'Optimizing path using {method} method')
    if method == 'IDPP':
        path = IDPPSolver.from_endpoints([is_pmg, fs_pmg], n_images, sort_tol=sort_tol)
        new_path = path.run(maxiter=5000, tol=1e-5, gtol=1e-3)
    elif method == 'linear':
        new_path = is_pmg.interpolate(fs_pmg, n_images+1, autosort_tol=sort_tol)
    else:
        raise ValueError(f'{method} not supported')
    
    # conver path to ase format, and add SinglePointCalculator
    ase_path = [i.to_ase_atoms() for i in new_path]
    ase_path[0].calc = SinglePointCalculator(ase_path[0].copy(), energy=is_e, forces=is_f)
    ase_path[-1].calc = SinglePointCalculator(ase_path[-1].copy(), energy=fs_e, forces=fs_f)

    # Set fix and magmom
    if bool(fix_height):
        for image in ase_path[1:-1]:
            set_fix_for_Atoms(image, fix_height, fix_dir)
    if bool(mag_ele):
        for image in ase_path[1:-1]:
            set_magmom_for_Atoms(image, mag_ele, mag_num)

    # write path
    print(f'Writing path: {output_file},Number of images: {len(ase_path)}')
    write(output_file, ase_path)

def main():

    # parse arguments
    parser = argparse.ArgumentParser(description='Make input files for NEB calculation')
    parser.add_argument('-n', type=int, help='Number of images', required=True)
    parser.add_argument('-f','--format', type=str, default='abacus-out', help='Format of the input files, default is abacus-out')
    parser.add_argument('-i', '--input',type=str, default=None, help='IS and FS file', nargs=2,required=True)
    parser.add_argument('-m', '--method',type=str, default='IDPP', help='Method to generate images', choices=['IDPP','linear'])
    parser.add_argument('-o', type=str, default='init_neb_chain.traj', help='Output file')
    parser.add_argument('-sort_tol', type=float, default=1.0, help='Sort tolerance for matching the initial and final structures, default is 1.0')
    parser.add_argument('--fix', type=str, default=None, help='[height]:[direction] : fix atom below height (fractional) in direction (0,1,2 for x,y,z), default None')
    parser.add_argument('--mag',type=str, default=None, help='[element1]:[magmom1],[element2]:[magmom2],... : set initial magmom for atoms of element, default None')

    args = parser.parse_args()

    fix_height = None
    fix_dir = None
    mag_ele = None
    mag_num = None

    if bool(args.fix):
        fix_height, fix_dir = args.fix.split(':')
        fix_height = float(fix_height)
        fix_dir = int(fix_dir)
    if bool(args.mag):
        mag_ele = [i.split(':')[0] for i in args.mag.split(',')]
        mag_num = [float(i.split(':')[1]) for i in args.mag.split(',')]
    
    generate(method=args.method, n_images=args.n, is_file=args.input[0], fs_file=args.input[1], 
             output_file=args.o, sort_tol=args.sort_tol, format=args.format,
             fix_height=fix_height, fix_dir=fix_dir,
             mag_ele=mag_ele, mag_num=mag_num)

if __name__ == '__main__':
    main()
