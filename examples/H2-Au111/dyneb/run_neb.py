# Run a example: H2 dissociation on Au(111) surface
# JamesMisaka in 2023-0925
# An example for read running_relax.log for init and final information
# Then do NEB calculation

import os 
from ase.calculators.abacus import AbacusProfile
from ase.optimize import FIRE, BFGS, QuasiNewton
from ase.io import read, write
from ase.constraints import FixAtoms
#from pathlib import Path

# set pythonpath: not useful
# ROOTPATH=os.path.abspath("../..")
# os.environ['PYTHONPATH'] = f'{ROOTPATH}:$PYTHONPATH'

from abacus_neb import AbacusNEB

# setting
directory = "OUT"
neb_optimizer = FIRE # suited for CI-NEB
algorism = "improvedtangent" # IT-NEB is recommended
#dyneb=True  # default
interpolate = "idpp" # linear or idpp
climb = True
n_max = 8
mpi = 4
omp = 8
abacus = 'abacus'
fix_height = 0.25 # define the height of fixed atoms
initial_result = 'H2/OUT.init/running_relax.log'
final_result = '2H/OUT.final/running_relax.log'
#example_dir="/lustre/home/2201110432/example/abacus"
example_dir="/data/home/liuzq/example/"
pseudo_dir = f"{example_dir}/PP"
basis_dir = f"{example_dir}/ORB"
pp = {
        "H": "H_ONCV_PBE-1.0.upf",
        "Au": "Au_ONCV_PBE-1.0.upf", 
        }
basis = {   
         "H": "H_gga_6au_100Ry_2s1p.orb",
        "Au": "Au_gga_7au_100Ry_4s2p2d1f.orb"
            }
kpts = [3, 3, 1]
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
    'efield_pos_max': 0.8,
}


# below is fixed workflow for do ASE-NEB-ABACUS

# set calculator
# os.environ['OMP_NUM_THREADS'] = f'{omp}'
# profile = AbacusProfile(
#     argv=['mpirun', '-np', f'{mpi}', abacus])

# Initial state read from ABACUS calculation result:
# reading problem is fixed in https://gitlab.com/1041176461/ase-abacus/-/commit/25e8f02cdfa23d70b2862417fb14457d9672a532
initial = read(initial_result, index=-1, format='abacus-out')

# Final state read frome ABACUS calculation result:
final = read(final_result, index=-1, format='abacus-out')

# should set fix in ASE itself, fix information cannot be read from abacus-out
mask = initial.get_scaled_positions()[:,2] < fix_height
fix = FixAtoms(mask=mask)
initial.set_constraint(fix)
final.set_constraint(fix)

# do neb calculation by DyNEB
neb = AbacusNEB(initial=initial, final=final, parameters=parameters,
                directory=directory, mpi=mpi, omp=omp, abacus=abacus, 
                algorism=algorism, n_max=n_max,)
neb.run(optimizer=neb_optimizer, climb=climb, interpolate=interpolate, fmax=0.05)

# Get barrier
barrier = neb.get_barriers()
print(barrier)
neb.plot_bands()

# Visualize the results
# os.system(f'ase gui neb.traj@-{n_max}:')

