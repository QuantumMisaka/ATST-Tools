# Run a example: Al diffusion on Cu(100) surface
# JamesMisaka in 2023-0919
# An example for read STRU and do init and final calculation
# Then do NEB calculation
# Example: N2 diffusion on Cu111 surface
# N2.stru and 2N.stru are generated from ASE

import os 
from ase.calculators.abacus import Abacus, AbacusProfile
from ase.optimize import FIRE, BFGS, QuasiNewton
from ase.io import read, write
#from pathlib import Path

# set pythonpath: not useful
# ROOTPATH=os.path.abspath("../..")
# os.environ['PYTHONPATH'] = f'{ROOTPATH}'

from abacus_neb import AbacusNEB

# setting
# optimizer = FIRE # suited for CI-NEB
init_directory = "INIT"
final_directory = "FINAL"
neb_directory = "OUT"
optimizer = FIRE
algorism = "improvedtangent" # IT-NEB is recommended
#dyneb=True  # default
interpolate = "linear" # linear or idpp
climb = True
n_max = 7
mpi = 1
omp = 16
abacus = "abacus"
#example_dir = "/lustre/home/2201110432/example"
#example_dir="/data/home/example/"
example_dir="/home/james/example"
pseudo_dir = f"{example_dir}/PP"
basis_dir = f"{example_dir}/ORB"
pp = {'N':'N_ONCV_PBE-1.0.upf',
      'Cu':'Cu_ONCV_PBE-1.0.upf',}
basis = {'N':'N_gga_7au_100Ry_2s2p1d.orb',
         'Cu':'Cu_gga_8au_100Ry_4s2p2d1f.orb',}
kpts = [3, 3, 1]
# for abacus to do relaxation
opt_parameters = {
    'calculation': 'relax',
    'xc': 'pbe',
    'ecutwfc': 100,
    'smearing_method': 'gaussian',
    'smearing_sigma': 0.002,
    'basis_type': 'lcao',
    'ks_solver': 'genelpa',
    'mixing_type': 'pulay',
    'scf_thr': 1e-6,
    'scf_nmax': 300,
    'relax_nmax': 400,
    'relax_method': "cg",
    'out_chg': 1,
    'out_bandgap': 1,
    'kpts': kpts,
    'pp': pp,
    'basis': basis,
    'pseudo_dir': pseudo_dir,
    'basis_dir': basis_dir,
    'vdw_method': 'd3_bj',
    'cal_force': 1,
    'cal_stress': 1,
    'out_stru': 1,
    'out_chg': 0,
    'out_bandgap': 0,
    'efield_flag': 1,
    'dip_cor_flag': 1,
    'efield_dir': 2,
    'efield_pos_max': 0.8,
}

neb_parameters = {
    'calculation': 'scf',
    'xc': 'pbe',
    'ecutwfc': 100,
    'smearing_method': 'gaussian',
    'smearing_sigma': 0.002,
    'basis_type': 'lcao',
    'ks_solver': 'genelpa',
    'mixing_type': 'pulay',
    'scf_thr': 1e-6,
    'scf_nmax': 300,
    'out_chg': 1,
    'out_bandgap': 1,
    'kpts': kpts,
    'pp': pp,
    'basis': basis,
    'pseudo_dir': pseudo_dir,
    'basis_dir': basis_dir,
    'vdw_method': 'd3_bj',
    'cal_force': 1,
    'cal_stress': 1,
    'out_stru': 1,
    'out_chg': 0,
    'out_bandgap': 0,
    'efield_flag': 1,
    'dip_cor_flag': 1,
    'efield_dir': 2,
    'efield_pos_max': 0.8,
}

os.environ['OMP_NUM_THREADS'] = f'{omp}'
profile = AbacusProfile(
    argv=['mpirun', '-np', f'{mpi}', abacus])

# Initial stru read from ABACUS, should do single point calculation
initial = read('./N2-Cu.stru', format='abacus')

# Final stru read frome ABACUS mshould do single point calculation
final = read('./2N-Cu.stru', format='abacus')

# do optimization
initial.calc = Abacus(profile=profile, directory=init_directory,
                    **opt_parameters)
initial.get_potential_energy()

final.calc = Abacus(profile=profile, directory=final_directory,
                    **opt_parameters)
final.get_potential_energy()

# do neb calculation by DyNEB
neb = AbacusNEB(initial=initial, final=final, parameters=neb_parameters,
                directory=neb_directory, mpi=mpi, omp=omp, abacus=abacus, 
                n_max=n_max,)
neb.run(optimizer=optimizer, climb=climb, interpolate=interpolate, fmax=0.05)

# Get barrier
barrier = neb.get_barriers()
print(barrier)
neb.plot_bands()

# Visualize the results
# os.system(f'ase gui neb.traj@-{n_max}:')

