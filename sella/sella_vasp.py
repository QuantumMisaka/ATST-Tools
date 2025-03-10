from ase.calculators.vasp import Vasp
from ase.io import read, write, Trajectory
from ase.visualize import view
from sella import Sella, Constraints
import os

ts = read("STRU", format="abacus")
pp_path = "/mnt/sg001/home/fz_pku_jh/PP_ORB/vasp/"
# abacus calculator setting
exec = 'vasp_std'
mpi=32
omp=1
kpts = [3, 1, 2]
task_type = "w"  # a for append and restart from trajectory, w for write new trajectory
parameters = {
    'encut': 520,
    'ispin': 2,
    'xc': 'pbe',
    'kpts': kpts,
    'isif': 2,
    'ismear': 1,
    'sigma': 0.1,
    'algo': 'fast',
    'amix': 0.2,
    'amix_mag': 0.6,
    'ncore': 8,
    'prec': 'accurate',
    'ediff': 1.0e-6,
    'nelm': 300,
    'nelmdl': -11,
}


# developers only
sella_eta = 0.005 # 0.002 ~ 0.01

if __name__ == "__main__":
    os.environ["VASP_PP_PATH"] = pp_path
    os.environ['OMP_NUM_THREADS'] = f'{omp}'

    ts.calc = Vasp(
        command=f"mpirun -np {mpi} {exec}",
        directory="VASP",
        **parameters,
        )

    # run from TS stru, but fmax more required
    # set cons is optional
    #cons = Constraints(ts_neb)
    #cons.fix_translation(ts_neb._get_constraints()[0].get_indices())
    traj = Trajectory("run_sella.traj", task_type)
    dyn = Sella(
        ts,
        trajectory=traj,
        eta = sella_eta,
    )
    dyn.run(fmax=0.05)

    # output TS stru file
    write("TS_sella.stru", ts, format="abacus")
    write("TS_sella.cif", ts, format="cif")