from ase.calculators.abacus import Abacus, AbacusProfile
from ase.io import read, write
from ase.visualize import view
from sella import Sella, Constraints
import os

ts = read("STRU", format="abacus")
# abacus calculator setting
abacus = "abacus"
mpi=64
omp=1
lib_dir = "/lustre/home/2201110432/example/abacus"
#lib_dir = "/data/home/liuzq/example"
pseudo_dir = f"{lib_dir}/PP"
basis_dir = f"{lib_dir}/ORB"
pp = {
        "H":"H_ONCV_PBE-1.0.upf",
        "C":"C_ONCV_PBE-1.0.upf",
        "O":"O_ONCV_PBE-1.0.upf",
        'Al': 'Al_ONCV_PBE-1.0.upf',
        'Si': 'Si_ONCV_PBE-1.0.upf',
        'Pd': 'Pd_ONCV_PBE-1.0.upf',
        'Au': 'Au_ONCV_PBE-1.0.upf',
      }
basis = {
        "H":"H_gga_6au_100Ry_2s1p.orb",
        "C":"C_gga_7au_100Ry_2s2p1d.orb",
        "O":"O_gga_7au_100Ry_2s2p1d.orb",
        'Al': 'Al_gga_7au_100Ry_4s4p1d.orb',
        'Si': 'Si_gga_7au_100Ry_2s2p1d.orb',
        'Pd': 'Pd_gga_7au_100Ry_4s2p2d1f.orb',
        'Au': 'Au_gga_7au_100Ry_4s2p2d1f.orb',
         }
kpts = [2, 2, 2]
parameters = {
    'calculation': 'scf',
    'nspin': 2,
    'xc': 'pbe',
    'ecutwfc': 100,
    'ks_solver': 'genelpa',
    'symmetry': 0,
    'vdw_method': 'd3_bj',
    'smearing_method': 'gau',
    'smearing_sigma': 0.002,
    'basis_type': 'lcao',
    'mixing_type': 'broyden',
    'mixing_ndim': 8,
    'scf_thr': 1e-7,
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
    'out_chg': -1,
    'out_mul': 1,
    'out_wfc_lcao': 0,
    'out_bandgap': 0,
}


# developers only
sella_eta = 0.005

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

    ts.calc = set_abacus_calc(abacus, parameters, f"ABACUS", mpi, omp)

    # run from TS stru, but fmax more required
    # set cons is optional
    #cons = Constraints(ts_neb)
    #cons.fix_translation(ts_neb._get_constraints()[0].get_indices())
    dyn = Sella(
        ts,
        trajectory='run_sella.traj',
        eta = sella_eta,
    )
    dyn.run(fmax=0.05)

    # output TS stru file
    write("TS_sella.stru", ts, format="abacus")
    write("TS_sella.cif", ts, format="cif")
