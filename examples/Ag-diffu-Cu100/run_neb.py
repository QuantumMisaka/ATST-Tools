# Run a example: Al diffusion on Cu(100) surface
# JamesMisaka in 2023-0920 update
# An example for read running_relax.log for init and final information
# Then do NEB calculation

import os 
from ase.calculators.abacus import AbacusProfile
from ase.optimize import FIRE, BFGS, QuasiNewton
from ase.io import read, write
#from pathlib import Path

# set pythonpath: not useful
# ROOTPATH=os.path.abspath("../..")
# os.environ['PYTHONPATH'] = f'{ROOTPATH}'

from abacus_neb import AbacusNEB

# setting
directory = 'OUT'
optimizer = FIRE # suited for CI-NEB
#optimizer = QuasiNewton
algorism = "improvedtangent" # IT-NEB is recommended
interpolate = "idpp" # linear or idpp
#dyneb=True  # default
climb = True
n_max = 7
mpi = 4
omp = 16
abacus = 'abacus'
pseudo_dir = "/data/home/liuzq/example/PP"
basis_dir = "/data/home/liuzq/example/ORB"
pp = {"Ag": "Ag_ONCV_PBE-1.0.upf",
        "Cu": "Cu_ONCV_PBE-1.0.upf", }
basis = {"Ag": "Ag_gga_7au_100Ry_4s2p2d1f.orb",
            "Cu": "Cu_gga_8au_100Ry_4s2p2d1f.orb"}
kpts = [4, 4, 1]
parameters = {
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
    'out_chg': 0,
    'out_bandgap': 0,
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
    'efield_flag': 0,
    'dip_cor_flag': 0,
    'efield_dir': 2,
    'efield_pos_max': 0.6,
}
os.environ['OMP_NUM_THREADS'] = f'{omp}'
profile = AbacusProfile(
    argv=['mpirun', '-np', f'{mpi}', abacus])

# Initial state read from ABACUS calculation result:
initial = read('init/OUT.A/running_relax.log', index=-1, format='abacus-out')

# Final state read frome ABACUS calculation result:
final = read('final/OUT.B/running_relax.log', index=-1, format='abacus-out')

# do neb calculation by DyNEB
neb = AbacusNEB(initial=initial, final=final, parameters=parameters,
                directory=directory, mpi=mpi, omp=omp, abacus=abacus, 
                algorism=algorism, n_max=n_max,)
neb.run(optimizer=optimizer, climb=climb, interpolate=interpolate, fmax=0.05)

# Get barrier
barrier = neb.get_barriers()
print(barrier)
neb.plot_bands()

# Visualize the results
# os.system(f'ase gui neb.traj@-{n_max}:')

