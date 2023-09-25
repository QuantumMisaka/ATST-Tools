# Run a example: Li diffusion on in Si
# JamesMisaka in 2023-0919
# An example for read STRU and do init and final calculation
# Then do NEB calculation
# Example from Colomb Academy for ABACUS catalysis practice

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
interpolate = "idpp" # linear or idpp
climb = True
n_max = 7
mpi = 1
omp = 16
abacus = 'abacus'
# example_dir = "/lustre/home/2201110432/example/abacus"
# pseudo_dir = f"{example_dir}/PP"
# basis_dir = f"{example_dir}/ORB"
pseudo_dir = "."
basis_dir = "."
pp = {"Li": "Li_ONCV_PBE-1.2.upf",
        "Si": "Si_ONCV_PBE-1.2.upf", }
basis = {"Li": "Li_gga_8au_100Ry_4s1p.orb",
            "Si": "Si_gga_8au_100Ry_2s2p1d.orb"}
kpts = [2, 2, 2]
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
    'kpts': kpts,
    'pp': pp,
    'basis': basis,
    'pseudo_dir': pseudo_dir,
    'basis_dir': basis_dir,
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

# Initial stru read from ABACUS, should do single point calculation
initial = read('./initial_stru', format='abacus')

# Final stru read frome ABACUS mshould do single point calculation
final = read('./final_stru', format='abacus')

# relax calculation by abacus
initial.calc = Abacus(profile=profile, directory=init_directory,
                    **parameters)
qn_init = optimizer(initial, trajectory='init_opt.traj')
qn_init.run(fmax=0.05)

final.calc = Abacus(profile=profile, directory=final_directory,
                    **parameters)
qn_final = optimizer(initial, trajectory='final_opt.traj')
qn_final.run(fmax=0.05)

# do neb calculation by DyNEB
neb = AbacusNEB(initial=initial, final=final, parameters=parameters,
                directory=neb_directory, mpi=mpi, omp=omp, abacus=abacus, 
                algorism=algorism, n_max=n_max, )
neb.run(optimizer=optimizer, climb=climb, 
        interpolate=interpolate, fmax=0.05)

# Get barrier
barrier = neb.get_barriers()
print(barrier)
neb.plot_bands()

# Visualize the results
# os.system(f'ase gui neb.traj@-{n_max}:')

