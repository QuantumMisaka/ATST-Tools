# JamesMisaka in 2023-11-05
# Run NEB calculation on NEB images by ASE-ABACUS
# part of ASE-NEB-ABACUS scripts

import os
from ase.calculators.abacus import Abacus, AbacusProfile
from ase.optimize import FIRE, BFGS
from ase.mep.neb import NEB, DyNEB # newest ase
from ase.io import read, write
from ase.parallel import world, parprint, paropen
from abacus_neb import AbacusNEB
#from pathlib import Path

# setting for NEB
mpi = 16
omp = 2
neb_optimizer = FIRE # suited for CI-NEB
neb_directory = "NEBrun"
algorism = "improvedtangent" # IT-NEB is recommended
neb_type = "neb"  # neb, dyneb, but not autoneb
init_chain = "init_neb_chain.traj"
climb = True
fmax = 0.05  # eV / Ang
parallel = True

# setting for calculator
abacus = "abacus"
lib_dir = "/lustre/home/2201110432/example/abacus"
#lib_dir = ""
pseudo_dir = f"{lib_dir}/PP"
basis_dir = f"{lib_dir}/ORB"
pp = {
        'H':'H_ONCV_PBE-1.0.upf',
        'C':'C_ONCV_PBE-1.0.upf',
        'O':'O_ONCV_PBE-1.0.upf',
        'Fe': 'Fe_ONCV_PBE-1.0.upf'
      }
basis = {
        'H':'H_gga_6au_100Ry_2s1p.orb',
        'C':'C_gga_7au_100Ry_2s2p1d.orb',
        'O':'O_gga_7au_100Ry_2s2p1d.orb',
        'Fe':'Fe_gga_8au_100Ry_4s2p2d1f.orb',
        }

kpts = [3, 1, 2]
parameters = {
    'calculation': 'scf',
    'nspin': 2,
    'xc': 'pbe',
    'ecutwfc': 100,
    'ks_solver': 'genelpa',
    'symmetry': 1,
    'symmetry_autoclose': 1,
    'vdw_method': 'none',
    'smearing_method': 'mp',
    'smearing_sigma': 0.002,
    'basis_type': 'lcao',
    'mixing_type': 'pulay',
    'scf_thr': 1e-6,
    'scf_nmax': 500,
    'kpts': kpts,
    'pp': pp,
    'basis': basis,
    'pseudo_dir': pseudo_dir,
    'basis_dir': basis_dir,
    'cal_force': 1,
    'cal_stress': 1,
    'out_stru': 1,
    'out_chg': 0,
    'out_bandgap': 1,
    'out_mul': 1,
    'efield_flag': 1,
    'dip_cor_flag': 1,
    'efield_dir': 1,
    'efield_pos_max': 0.7
}


if __name__ == "__main__":
# running process
# read initial guessed neb chain
    init_chain = read(init_chain, index=':')
    neb = AbacusNEB(init_chain, parameters=parameters, parallel=parallel,
                    directory=neb_directory, mpi=mpi, omp=omp, abacus=abacus, 
                    algorism=algorism)
    neb.run(optimizer=neb_optimizer, climb=climb, fmax=fmax)


