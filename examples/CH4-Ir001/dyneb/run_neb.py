# Run a example: CH4 dissociation on Ir(001) surface
# JamesMisaka in 2023-0925, ref by Sob 1KFP NEB6
# An example for read calculation result by ABACUS
# Then do NEB calculation by ASE-ABACUS


import os 
from ase.calculators.abacus import Abacus, AbacusProfile
from ase.optimize import FIRE, BFGS, QuasiNewton
from ase.io import read, write
from ase.constraints import FixAtoms

# set pythonpath: not useful
# ROOTPATH=os.path.abspath("../..")
# os.environ['PYTHONPATH'] = f'{ROOTPATH}'

from abacus_neb import AbacusNEB

# setting
neb_optimizer = FIRE # suited for CI-NEB
init_directory = "INIT"
final_directory = "FINAL"
neb_directory = "OUT"
algorism = "improvedtangent" # IT-NEB is recommended
#dyneb=True  # default
interpolate = "idpp" # linear or idpp
climb = True
n_max = 8
mpi = 1
omp = 32
abacus = "abacus"
fix_height = 0.20 # define the height of fixed atoms
initial_result = 'init/OUT.init/running_relax.log'
final_result = 'final/OUT.final/running_relax.log'
#example_dir = "/lustre/home/2201110432/example/abacus"
example_dir="/data/home/liuzq/example/"
#example_dir="/home/james/example"
pseudo_dir = f"{example_dir}/PP"
basis_dir = f"{example_dir}/ORB"
pp = {
    'Ir':'Au_ONCV_PBE-1.0.upf',
      'H':'H_ONCV_PBE-1.0.upf',
      'C':'C_ONCV_PBE-1.0.upf',
      'O': 'O_ONCV_PBE-1.0.upf',
      }
basis = {
    'Ir':'Ir_gga_7au_100Ry_4s2p2d1f.orb',
         'H':'H_gga_6au_100Ry_2s1p.orb',
         'C': 'C_gga_7au_100Ry_2s2p1d.orb',
         'O': 'O_gga_7au_100Ry_2s2p1d.orb'
         ,}
kpts = [4, 4, 1]
# for abacus to do relaxation

parameters = {
    'calculation': 'scf',
    'xc': 'pbe',
    'ecutwfc': 100,
    'smearing_method': 'mp',
    'smearing_sigma': 0.008,
    'basis_type': 'lcao',
    'ks_solver': 'genelpa',
    'mixing_type': 'pulay',
    'scf_thr': 1e-6,
    'scf_nmax': 200,
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
    'efield_pos_max': 0.7
}

# below is fixed workflow for do ASE-NEB-ABACUS

# set calculator

os.environ['OMP_NUM_THREADS'] = f'{omp}'
profile = AbacusProfile(
    argv=['mpirun', '-np', f'{mpi}', abacus])

# Initial stru read from ABACUS
# reading problem is fixed in https://gitlab.com/1041176461/ase-abacus/-/commit/25e8f02cdfa23d70b2862417fb14457d9672a532
initial = read(initial_result,index=-1, format='abacus-out')

# Final stru read from ABACUS
final = read(final_result,index=-1, format='abacus-out')

# should set fix in ASE itself, fix information cannot be read from abacus-out
fix_indices = [atom.index for atom in initial
               [initial.get_scaled_positions()[:,2] < fix_height]]
fix = FixAtoms(indices=fix_indices)
initial.set_constraint(fix)
final.set_constraint(fix)

# do neb calculation by DyNEB
neb = AbacusNEB(initial=initial, final=final, parameters=parameters,
                directory=neb_directory, mpi=mpi, omp=omp, abacus=abacus, 
                algorism=algorism, n_max=n_max,)
neb.run(optimizer=neb_optimizer, climb=climb, interpolate=interpolate, fmax=0.05)

# Get barrier
barrier = neb.get_barriers()
print(barrier)
neb.plot_bands()

# Visualize the results
# os.system(f'ase gui neb.traj@-{n_max}:')

