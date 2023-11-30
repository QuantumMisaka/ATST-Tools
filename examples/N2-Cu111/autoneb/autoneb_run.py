# JamesMisaka in 2023-11-05
# Run AutoNEB calculation on NEB images by ASE-ABACUS
# part of ASE-NEB-ABACUS scripts

from ase.optimize import FIRE, BFGS
from ase.io import read, write
from ase.parallel import world, parprint, paropen
#from pathlib import Path
from abacus_autoneb import AbacusAutoNEB

# setting for NEB
mpi = 16
omp = 4
neb_optimizer = FIRE # suited for CI-NEB
neb_directory = "AutoNEBrun"
# algorism = 'eb' # default for AutoNEB
algorism = "improvedtangent" # IT-NEB
init_chain = "init_neb_chain.traj"
climb = True
fmax = [0.20, 0.05]  # eV / Ang, 2 fmax for others and last CI-NEB in all-images
n_simul = world.size # only for autoneb, number of simultaneous calculation
n_images = 10 # only for autoneb, max number of all image, which should be reached

# setting for calculator
abacus = "abacus"
lib_dir = "/lustre/home/2201110432/example/abacus"
#lib_dir = ""
pseudo_dir = f"{lib_dir}/PP"
basis_dir = f"{lib_dir}/ORB"
pp = {
        'N':'N_ONCV_PBE-1.0.upf',
        'Cu':'Cu_ONCV_PBE-1.0.upf',
      }
basis = {
        'N':'N_gga_7au_100Ry_2s2p1d.orb',
        'Cu':'Cu_gga_8au_100Ry_4s2p2d1f.orb',
        }
kpts = [3, 1, 3]
parameters = {
    'calculation': 'scf',
    'nspin': 2,
    'xc': 'pbe',
    'ecutwfc': 100,
    'ks_solver': 'genelpa',
    'symmetry': 0,
    'vdw_method': 'd3_bj',
    'smearing_method': 'gaussian',
    'smearing_sigma': 0.001,
    'basis_type': 'lcao',
    'mixing_type': 'broyden',
    'scf_thr': 1e-6,
    'scf_nmax': 100,
    'kpts': kpts,
    'pp': pp,
    'basis': basis,
    'pseudo_dir': pseudo_dir,
    'basis_dir': basis_dir,
    'cal_force': 1,
    'cal_stress': 1,
    'out_stru': 1,
    'out_chg': 0,
    'out_mul': 0,
    'out_bandgap': 0,
    'efield_flag': 1,
    'dip_cor_flag': 1,
    'efield_dir': 1,
    'efield_pos_max': 0.0
}


if __name__ == "__main__": 
# running process
# read initial guessed neb chain
    init_chain = read(init_chain, index=':')
    neb = AbacusAutoNEB(init_chain, parameters, algorism=algorism, 
                        directory=neb_directory,
                        n_simul=n_simul, n_max=n_images, 
                        abacus=abacus,  mpi=mpi, omp=omp, )
    neb.run(optimizer=neb_optimizer, climb=climb, fmax=fmax)
