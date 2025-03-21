from ase.calculators.abacus import Abacus, AbacusProfile
from ase.io import read, write, Trajectory
from ase.visualize import view
from sella import IRC
import os

ts_opt = read("STRU", format="abacus")
# abacus calculator setting
abacus = "abacus"
mpi=16
omp=4
lib_dir = "/lustre/home/2201110432/example/abacus"
task_type = "w"  # a for append and restart from trajectory, w for write new trajectory
#lib_dir = "/data/home/liuzq/example"
pseudo_dir = f"{lib_dir}/PP"
basis_dir = f"{lib_dir}/ORB"
pp = {
        "H":"H_ONCV_PBE-1.0.upf",
        "C":"C_ONCV_PBE-1.0.upf",
        "O":"O_ONCV_PBE-1.0.upf",
        "Fe":"Fe_ONCV_PBE-1.0.upf",
      }
basis = {
        "H":"H_gga_6au_100Ry_2s1p.orb",
        "C":"C_gga_7au_100Ry_2s2p1d.orb",
        "O":"O_gga_7au_100Ry_2s2p1d.orb",
        "Fe":"Fe_gga_8au_100Ry_4s2p2d1f.orb",
         }
kpts = [3, 1, 2]
irc_log = "irc_log.traj"
dx = 0.1
fmax = 0.05
steps = 1000
parameters = {
    'calculation': 'scf',
    'nspin': 2,
    'xc': 'pbe',
    'ecutwfc': 100,
    'ks_solver': 'genelpa',
    'symmetry': 0,
    'vdw_method': 'none',
    'smearing_method': 'mp',
    'smearing_sigma': 0.002,
    'basis_type': 'lcao',
    'mixing_type': 'broyden',
    'mixing_ndim': 20,
    'scf_thr': 1e-7,
    'scf_nmax': 300,
    'kpts': kpts,
    'pp': pp,
    'basis': basis,
    'pseudo_dir': pseudo_dir,
    'basis_dir': basis_dir,
    'cal_force': 1,
    'cal_stress': 0,
    'init_wfc': 'atomic',
    'init_chg': 'atomic',
    'out_stru': 1,
    'out_chg': -1,
    'out_mul': 1,
    'out_wfc_lcao': 0,
    'out_bandgap': 0,
    'efield_flag': 1,
    'dip_cor_flag': 1,
    'efield_dir': 1,
}

# developers only
sella_eta = 0.002


# set calculator
def set_abacus_calc(abacus, parameters, directory, mpi, omp) -> Abacus:
    """Set Abacus calculators"""
    os.environ['OMP_NUM_THREADS'] = f'{omp}'
    profile = AbacusProfile(f"mpirun -np {mpi} {abacus}")
    out_directory = directory
    calc = Abacus(profile=profile, directory=out_directory,
                **parameters)
    return calc

if __name__ == "__main__":
    ts_opt.calc = set_abacus_calc(abacus, parameters, f"ABACUS", mpi, omp)
    # run from TS stru, but fmax more required
    # set cons is optional
    #cons = Constraints(ts_neb)
    #cons.fix_translation(ts_neb._get_constraints()[0].get_indices())
    irc_traj = Trajectory(irc_log, task_type)
    irc = IRC(ts_opt, trajectory=irc_traj, dx=dx, eta=sella_eta)
    irc.run(fmax, steps=steps, direction='forward')
    irc.run(fmax, steps=steps, direction='reverse')
    # normalize the trajectory
    irc_log_norm = []
    irc_origin = read(irc_log,":")
    ene_last = -99999999
    for stru in irc_origin[::-1]:
        if stru.get_potential_energy() > ene_last:
            irc_log_norm.append(stru)
            ene_last = stru.get_potential_energy()
        else:
            break
    ene_last = 99999999
    for stru in irc_origin:
        if stru.get_potential_energy() < ene_last:
            irc_log_norm.append(stru)
            ene_last = stru.get_potential_energy()
        else:
            break
    write(f"norm_{irc_log}", irc_log_norm, format="traj")