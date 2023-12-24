# read STRU file, generate INCAR and KPT file 
# run optimization jobs
import os
from ase.optimize import QuasiNewton, BFGS, LBFGS, FIRE, GPMin, MDMin
from ase.io import read
from ase.io import Trajectory
from ase.calculators.abacus import Abacus, AbacusProfile

# setting
stru = read('STRU', format='abacus')
optimizer = BFGS
mpi = 16
omp = 4
abacus = 'abacus'
lib_dir = "/lustre/home/2201110432/example/abacus"
#lib_dir = "/data/home/liuzq/example"
pseudo_dir = f"{lib_dir}/PP"
basis_dir = f"{lib_dir}/ORB"
out_dir = "OUT"
properties = ["energy", "forces", "stress"]
out_traj = Trajectory("relax.traj", 'w', stru, properties=properties)
pp = {
        'H':'H_ONCV_PBE-1.0.upf',
        'C':'C_ONCV_PBE-1.0.upf',
        'O':'O_ONCV_PBE-1.0.upf',
        'Fe':'Fe_ONCV_PBE-1.0.upf',
      }
basis = {
        'H':'H_gga_6au_100Ry_2s1p.orb',
        'C': 'C_gga_7au_100Ry_2s2p1d.orb',
        'O': 'O_gga_7au_100Ry_2s2p1d.orb',
        'Fe': 'Fe_gga_8au_100Ry_4s2p2d1f.orb',
        }
kpts = [3, 1, 2] # KPT setting (will be generate next)
# INPUT setting
parameters = {
    'calculation': 'scf',
    'basis_type': 'lcao',
    'ks_solver': 'genelpa',
    'vdw_method': 'none',
    'nspin' : 2,
    'xc': 'pbe',
    'ecutwfc': 100,
    'kpts': kpts,
    'pp': pp,
    'basis': basis,
    'pseudo_dir': pseudo_dir,
    'basis_dir': basis_dir,
    'smearing_method': 'mp',
    'smearing_sigma': 0.002,
    'mixing_type': 'broyden',
    'mixing_beta': 0.4,
    'mixing_gg0': 1.0,
    'mixing_ndim': 8,
    'scf_thr': 1e-6,
    'init_wfc': 'file',
    'cal_force': 1,
    'cal_stress': 1,
    'out_stru': 1,
    'out_chg': 0,
    'out_bandgap': 1,
    'out_wfc_lcao': 1,
    'efield_flag': 1,
    'dip_cor_flag': 1,
    'efield_dir': 1,
    'efield_pos_max': 0.7,
}

# running
os.environ['OMP_NUM_THREADS'] = f'{omp}'
profile = AbacusProfile(
    argv=['mpirun', '-np', f'{mpi}', abacus])
stru.calc = Abacus(profile=profile, directory=out_dir,
                    **parameters)
qn = optimizer(stru, trajectory=out_traj)
qn.run(fmax=0.05)
