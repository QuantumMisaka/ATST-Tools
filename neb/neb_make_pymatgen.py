import argparse
from ase.io import read, write
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.diffusion.neb.pathfinder import IDPPSolver
from ase.calculators.singlepoint import SinglePointCalculator

'''
Author: mosey
Last update: 2024-03-28

Dependencies:
1. Install pymatgen: pip install pymatgen
2. Install pymatgen-analysis-diffusion: pip install pymatgen-analysis-diffusion
3. Install ase: pip install ase
'''

def generate(method:str, n_images:int, is_file:str, fs_file:str, 
             output_file:str,sort_tol:float,format:str):

    # read files
    print(f'Reading files: {is_file} and {fs_file}')
    is_atom = read(is_file,format=format)
    fs_atom = read(fs_file,format=format)
    is_e = is_atom.get_potential_energy()
    fs_e = fs_atom.get_potential_energy()
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
    ase_path[0].calc = SinglePointCalculator(ase_path[0].copy(), energy=is_e)
    ase_path[-1].calc = SinglePointCalculator(ase_path[-1].copy(), energy=fs_e)
    print(f'Writing path: {output_file},Number of images: {len(ase_path)}')
    write(output_file, ase_path)

def main():
    parser = argparse.ArgumentParser(description='Make input files for NEB calculation')
    parser.add_argument('-n', type=int, help='Number of images', required=True)
    parser.add_argument('-f','--format', type=str, default='abacus-out', help='Format of the input files, default is abacus-out')
    parser.add_argument('-i', '--input',type=str, default=None, help='IS and FS file', nargs=2,required=True)
    parser.add_argument('-m', '--method',type=str, default='IDPP', help='Method to generate images', choices=['IDPP','linear'])
    parser.add_argument('-o', type=str, default='init_neb_chain.traj', help='Output file')
    parser.add_argument('-sort_tol', type=float, default=1.0, help='Sort tolerance for matching the initial and final structures, default is 1.0')
    args = parser.parse_args()

    generate(args.method, args.n, args.input[0], args.input[1], args.o, args.sort_tol, args.format)

if __name__ == '__main__':
    main()