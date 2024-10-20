from ase.optimize import FIRE, BFGS
from ase.io import read, write
from abacus_neb import AbacusNEB

# setting
mpi = 16
omp = 1
fmax = 0.05  # eV / Ang
neb_optimizer = FIRE # suited for CI-NEB
neb_directory = "NEBrun"
algorism = "improvedtangent" # IT-NEB is recommended
climb = True
dyneb = True
parallel = False
k = 0.10 # eV/Ang^2, spring constant
init_chain = "init_neb_chain.traj"
abacus = "abacus"
#lib_dir = "/lustre/home/2201110432/example/abacus"
lib_dir = ""
pseudo_dir = f"{lib_dir}"
basis_dir = f"{lib_dir}"
pp = {
      'Li':'C_ONCV_PBE-1.2.upf',
      'Si':'H_ONCV_PBE-1.2.upf',
      }
basis = {
         'Li': 'Li_gga_8au_100Ry_4s1p.orb',
         'Si': 'Si_gga_8au_100Ry_2s2p1d.orb',
         ,}
kpts = [2, 2, 2]
parameters = {
    'calculation': 'scf',
    'nspin': 1,
    'xc': 'pbe',
    'ecutwfc': 100,
    'dft_functional': 'pbe',
    'ks_solver': 'genelpa',
    'symmetry': 0,
    'vdw_method': 'none',
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
    'init_wfc': 'atomic',
    'init_chg': 'atomic',
    'cal_force': 1,
    'cal_stress': 1,
    'out_stru': 1,
    'out_chg': -1,
}


if __name__ == "__main__":
# running process
# read initial guessed neb chain
    init_chain = read(init_chain, index=':')
    neb = AbacusNEB(init_chain, parameters=parameters, parallel=parallel,
                    directory=neb_directory, mpi=mpi, omp=omp, abacus=abacus, 
                    algorism=algorism, k=k, dyneb=dyneb)