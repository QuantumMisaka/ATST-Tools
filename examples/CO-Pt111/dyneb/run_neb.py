# Run a example: CO dissociation on Pt(111) surface
# JamesMisaka in 2023-0924
# An example for read STRU to do init and final optimization with ASE
# Then do NEB calculation with ASE

import os 
from ase.calculators.abacus import Abacus, AbacusProfile
from ase.optimize import FIRE, BFGS, QuasiNewton
from ase.io import read, write
#from ase.constraints import FixAtoms
#from pathlib import Path

# set pythonpath: not useful
# ROOTPATH=os.path.abspath("../..")
# os.environ['PYTHONPATH'] = f'{ROOTPATH}:$PYTHONPATH'

from abacus_neb import AbacusNEB

# setting
neb_optimizer = FIRE # suited for CI-NEB
opt_optimizer = QuasiNewton
algorism = "improvedtangent" # IT-NEB is recommended
#dyneb=True  # default
interpolate = "idpp" # linear or idpp
climb = True
n_max = 4
mpi = 1
omp = 16
abacus = 'abacus'
init_stru = 'init/STRU'
final_stru = 'final/STRU'
neb_directory = "OUT"
init_directory = "INIT"
final_directory = "FINAL"
#example_dir="/lustre/home/2201110432/example/abacus"
example_dir="/data/home/liuzq/example/"
pseudo_dir = f"{example_dir}/PP"
basis_dir = f"{example_dir}/ORB"
pp = {
        "C": "C_ONCV_PBE-1.0.upf",
        "O": "O_ONCV_PBE-1.0.upf",
        "Pt": "Pt_ONCV_PBE-1.0.upf", 
        }
basis = {
        "C": "C_gga_7au_100Ry_2s2p1d.orb",
        "O": "O_gga_7au_100Ry_2s2p1d.orb",
        "Pt": "Pt_gga_7au_100Ry_4s2p2d1f.orb"
            }
kpts = [3, 3, 1]
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
    'efield_pos_max': 0.7,
}
os.environ['OMP_NUM_THREADS'] = f'{omp}'
profile = AbacusProfile(
    argv=['mpirun', '-np', f'{mpi}', abacus])

# Initial state read 
initial = read(init_stru, format='abacus')

# relax calculation by ase-abacus, abacus do scf
initial.calc = Abacus(profile=profile, directory=init_directory,
                    **parameters)
qn_init = opt_optimizer(initial, trajectory='init_opt.traj')
qn_init.run(fmax=0.05)


# Final state read 
final = read(final_stru, format='abacus')
final.calc = Abacus(profile=profile, directory=final_directory,
                    **parameters)
qn_final = opt_optimizer(final, trajectory='final_opt.traj')
qn_final.run(fmax=0.05)

# fix is done by reading from STRU and do relax by ASE-ABACUS
# in OUT directory, fix is done in print-out STRU
# but it is in question whether the fix is effective ? need more test

# do neb calculation by DyNEB
neb = AbacusNEB(initial=initial, final=final, parameters=parameters,
                directory=neb_directory, mpi=mpi, omp=omp, abacus=abacus, 
                algorism=algorism, n_max=n_max,)
neb.run(optimizer=neb_optimizer, climb=climb, 
        interpolate=interpolate, fmax=0.05)

# Get barrier
barrier = neb.get_barriers()
print(barrier)
neb.plot_bands()

# Visualize the results
# os.system(f'ase gui neb.traj@-{n_max}:')

