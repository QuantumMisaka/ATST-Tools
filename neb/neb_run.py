# JamesMisaka in 2023-11-27
# Run NEB calculation on NEB images by ASE-ABACUS
# part of ATST-Tools scripts

from ase.optimize import FIRE, BFGS
from ase.io import read, write
from abacus_neb import AbacusNEB

# setting
mpi = 8
omp = 4
fmax = 0.05  # eV / Ang
neb_optimizer = FIRE # suited for CI-NEB
neb_directory = "NEBrun"
algorism = "improvedtangent" # IT-NEB is recommended
climb = True
dyneb = False  
parallel = True
k = 0.10 # eV/Ang^2, spring constant
init_chain = "init_neb_chain.traj"
abacus = "abacus"
#lib_dir = "/lustre/home/2201110432/example/abacus"
lib_dir = ""
pseudo_dir = f"{lib_dir}"
basis_dir = f"{lib_dir}"
pp = {
      'C':'C_ONCV_PBE-1.0.upf',
      'H':'H_ONCV_PBE-1.0.upf',
      'Pt':'Pt_ONCV_PBE-1.0.upf',
      }
basis = {
         'C': 'C_gga_7au_100Ry_2s2p1d.orb',
         'H': 'H_gga_6au_100Ry_2s1p.orb',
         'Pt': 'Pt_gga_7au_100Ry_4s2p2d1f.orb'
         ,}
kpts = [2, 1, 2]
parameters = {
    'calculation': 'scf',
    'nspin': 2,
    'xc': 'pbe',
    'ecutwfc': 100,
    'dft_functional': 'pbe',
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
    'init_wfc': 'atomic',
    'init_chg': 'atomic',
    'cal_force': 1,
    'cal_stress': 1,
    'out_stru': 1,
    'out_chg': 0,
    'out_mul': 0,
    'out_wfc_lcao': 0,
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
    neb = AbacusNEB(init_chain, parameters=parameters, parallel=parallel,
                    directory=neb_directory, mpi=mpi, omp=omp, abacus=abacus, 
                    algorism=algorism, k=k, dyneb=dyneb)
    neb.run(optimizer=neb_optimizer, climb=climb, fmax=fmax)


