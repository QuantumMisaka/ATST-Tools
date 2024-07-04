from ase.calculators.abacus import Abacus, AbacusProfile
from ase.io import read, write
from ase.visualize import view
from sella import Sella, Constraints
import os

# abacus calculator setting
abacus = "abacus"
#lib_dir = "/lustre/home/2201110432/example/abacus"
lib_dir = "/home/james/example/abacus"
pseudo_dir = f"{lib_dir}/PP"
basis_dir = f"{lib_dir}/ORB"
pp = {
      'Au':'Au_ONCV_PBE-1.0.upf',
      'H':'H_ONCV_PBE-1.0.upf',
      }
basis = {
         'Au': 'Au_gga_7au_100Ry_4s2p2d1f.orb',
         'H': 'H_gga_6au_100Ry_2s1p.orb',
         }
kpts = [3, 1, 3]
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
    'cal_force': 1,
    'cal_stress': 1,
    'init_wfc': 'atomic',
    'init_chg': 'atomic',
    'out_stru': 1,
    'out_chg': 0,
    'out_mul': 0,
    'out_wfc_lcao': 0,
    'out_bandgap': 0,
    'efield_flag': 1,
    'dip_cor_flag': 1,
    'efield_dir': 1,
}


# read traj of NEB
neb_traj = read("neb_H2Au111.traj", ":")
ts_neb = neb_traj[6]
ts_bef_neb = neb_traj[5]
ts_aft_neb = neb_traj[7]

# set calculator
def set_abacus_calc(abacus, parameters, directory, mpi, omp) -> Abacus:
    """Set Abacus calculators"""
    os.environ['OMP_NUM_THREADS'] = f'{omp}'
    profile = AbacusProfile(
        argv=['mpirun', '-np', f'{mpi}', abacus])
    out_directory = directory
    calc = Abacus(profile=profile, directory=out_directory,
                **parameters)
    return calc

ts_neb.calc = set_abacus_calc(abacus, parameters, "Sella-TS", 4, 4)
ts_bef_neb.calc = set_abacus_calc(abacus, parameters, "Sella-beforeTS", 4, 4)
ts_aft_neb.calc = set_abacus_calc(abacus, parameters, "Sella-afterTS", 4, 4)

# run from TS stru, but fmax more required
# set cons is optional
cons = Constraints(ts_neb)
cons.fix_translation(ts_neb._get_constraints()[0].get_indices())
print("==> run from TS stru <==")
dyn = Sella(
    ts_neb,
    constraints=cons,
    trajectory='H2-Au111-Sella-from_NEB_TS.traj',
)
dyn.run(fmax=0.03)

# run from TS_bef_neb with NEB-TS fmax
cons = Constraints(ts_bef_neb)
cons.fix_translation(ts_bef_neb._get_constraints()[0].get_indices()) # which seems to have problem
print("==> run from TS_before <==")
dyn = Sella(
    ts_neb,
    constraints=cons,
    trajectory='H2-Au111-Sella-from_NEB_TS_before.traj',
)
dyn.run(fmax=0.05)

# run from TS_aft_neb with NEB-TS fmax
cons = Constraints(ts_aft_neb)
cons.fix_translation(ts_aft_neb._get_constraints()[0].get_indices())
print("==> run from TS_after <==")
dyn = Sella(
    ts_neb,
    constraints=cons,
    trajectory='H2-Au111-Sella-from_NEB_TS_after.traj',
)
dyn.run(fmax=0.05)



