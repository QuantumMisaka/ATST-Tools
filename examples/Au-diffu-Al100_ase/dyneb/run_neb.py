# Run a example: Au diffusion on Al(100) surface
# JamesMisaka in 2023-0918

import os 
from ase.calculators.abacus import Abacus, AbacusProfile
from ase.optimize import FIRE, BFGS
from ase.io import read, write

# set pythonpath: not useful
# ROOTPATH=os.path.abspath("../..")
# os.environ['PYTHONPATH'] = f'{ROOTPATH}:$PYTHONPATH'

from abacus_neb import AbacusNEB

directory = 'OUT'
optimizer = BFGS # suited for CI-NEB
interpolate = "linear" # linear or idpp
n_max = 5
mpi = 64
omp = 1
abacus = 'abacus'
pseudo_dir = "/data/home/liuzq/example/PP"
basis_dir = "/data/home/liuzq/example/ORB"
pp = {"Al": "Al_ONCV_PBE-1.0.upf",
        "Au": "Au_ONCV_PBE-1.0.upf", }
basis = {"Al": "Al_gga_7au_100Ry_4s4p1d.orb",
            "Au": "Au_gga_7au_100Ry_4s2p2d1f.orb"}
kpts = [2, 2, 1]
parameters = {
    'calculation': 'scf',
    'xc': 'pbe',
    'ecutwfc': 100,
    'smearing_method': 'gaussian',
    'smearing_sigma': 0.01,
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
    'efield_pos_max': 0.6,
}

from ase.build import fcc100, add_adsorbate
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
slab = fcc100('Al', size=(2, 2, 3))
add_adsorbate(slab, 'Au', 1.7, 'hollow')
slab.center(axis=2, vacuum=4.0)
mask = [atom.tag > 1 for atom in slab]
slab.set_constraint(FixAtoms(mask=mask))
os.environ['OMP_NUM_THREADS'] = f'{omp}'
profile = AbacusProfile(
    argv=['mpirun', '-np', f'{mpi}', abacus])
slab.calc = Abacus(profile=profile, directory="INIT",
                    **parameters)

qn = QuasiNewton(slab, trajectory='initial.traj')
qn.run(fmax=0.05)
slab[-1].x += slab.get_cell()[0, 0] / 2
qn = QuasiNewton(slab, trajectory='final.traj')
qn.run(fmax=0.05)
# Initial state:
initial = read('initial.traj')

# Final state:
final = read('final.traj')

neb = AbacusNEB(initial=initial, final=final, parameters=parameters,
                directory=directory,
                mpi=mpi, omp=omp, abacus=abacus, n_max=n_max)
neb.run(optimizer=optimizer, climb=False, interpolate=interpolate, fmax=0.05)

# Get barrier
barrier = neb.get_barriers()
print(barrier)
neb.plot_bands()

# Visualize the results
# os.system(f'ase gui neb.traj@-{n_max}:')

