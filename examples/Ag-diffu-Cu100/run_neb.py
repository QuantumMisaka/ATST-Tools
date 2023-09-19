# Run a example: Al diffusion on Cu(100) surface
# JamesMisaka in 2023-0919
# Should do ABACUS calculation for init and final state first

import os 
from ase.calculators.abacus import AbacusProfile
from ase.optimize import FIRE, BFGS
from ase.io import read, write
from ase.optimize import QuasiNewton
#from pathlib import Path

# set pythonpath
ROOTPATH=os.path.abspath("../..")
os.environ['PYTHONPATH'] = f'{ROOTPATH}'

from abacus_neb import AbacusNEB

# setting
n_max = 7
mpi = 4
omp = 16
abacus = 'abacus'
pseudo_dir = "/data/home/liuzq/example/PP"
basis_dir = "/data/home/liuzq/example/ORB"
pp = {"Ag": "Al_ONCV_PBE-1.0.upf",
        "Cu": "Cu_ONCV_PBE-1.0.upf", }
basis = {"Ag": "Ag_gga_7au_100Ry_4s4p1d.orb",
            "Cu": "Cu_gga_7au_100Ry_4s2p2d1f.orb"}
kpts = [4, 4, 1]
parameters = {
    'calculation': 'scf',
    'xc': 'pbe',
    'ecutwfc': 100,
    'smearing_method': 'gaussian',
    'smearing_sigma': 0.004,
    'basis_type': 'lcao',
    'ks_solver': 'genelpa',
    'mixing_type': 'pulay',
    'scf_thr': 1e-6,
    'out_chg': 1,
    'out_bandgap': 1,
    'kpts': kpts,
    'pp': pp,
    'basis': basis,
    'pseudo_dir': pseudo_dir,
    'basis_dir': basis_dir,
    'vdw_method': 'd3_bj',
    'cal_force': 1,
    
}
os.environ['OMP_NUM_THREADS'] = f'{omp}'
profile = AbacusProfile(
    argv=['mpirun', '-np', f'{mpi}', abacus])

# Initial state read from ABACUS calculation result:
initial = read('OUT.init/running_relax.log', index=-1, format='abacus-out')

# Final state read frome ABACUS calculation result:
final = read('OUT.final/running_relax.log', index=-1, format='abacus-out')

# do neb calculation by DyNEB
neb = AbacusNEB(initial=initial, final=final, parameters=parameters,
                mpi=mpi, omp=omp, abacus=abacus, n_max=n_max)
neb.run(optimizer=BFGS, climb=False, interpolate=True, fmax=0.05)

# Get barrier
barrier = neb.get_barriers()
print(barrier)
neb.plot_bands()

# Visualize the results
# os.system(f'ase gui neb.traj@-{n_max}:')

