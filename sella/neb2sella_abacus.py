# neb+sella TS-Search by ABACUS
import numpy as np
import os, sys

from ase.io import read, write, Trajectory
from ase import Atoms
from ase.optimize import BFGS, FIRE, QuasiNewton
from ase.mep.neb import NEBTools, NEB, DyNEB
from ase.mep.autoneb import AutoNEB
from sella import Sella, Constraints
from ase.calculators.abacus import Abacus, AbacusProfile
from copy import deepcopy

from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.diffusion.neb.pathfinder import IDPPSolver


n_neb_images = 8
neb_fmax = 0.80  # neb should be rough
sella_fmax = 0.05 # dimer use neb guess
climb = True
neb_algorism = "improvedtangent"
neb_log = "neb_images.traj"
sella_log = "sella_images.traj"
OPTSolver = QuasiNewton
NEBSolver = FIRE

# abacus setting
abacus = "abacus"
mpi = 16
omp = 4
abacus = "abacus"
lib_dir = "/lustre/home/2201110432/example/abacus"
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
# INPUT setting
parameters = {
    'calculation': 'scf',
    'basis_type': 'lcao',
    'xc': 'pbe',
    'ks_solver': 'genelpa',
    'vdw_method': 'none',
    'nspin' : 2,
    'ecutwfc': 100,
    'pp': pp,
    'basis': basis,
    'pseudo_dir': pseudo_dir,
    'basis_dir': basis_dir,
    'smearing_method': 'mp',
    'smearing_sigma': 0.002,
    'mixing_type': 'broyden',
    'mixing_beta': 0.4,
    'mixing_gg0': 1.0,
    'mixing_ndim': 20,
    'scf_thr': 1e-7,
    'scf_nmax': 300,
    'kpts': kpts,
    'init_wfc': 'atomic',
    'init_chg': 'atomic',
    'cal_force': 1,
    'cal_stress': 0,
    'out_stru': 1,
    'out_chg': -1,
    'out_mul': 1,
    'out_bandgap': 1,
    'out_wfc_lcao': 0,
    'efield_flag': 1,
    'dip_cor_flag': 1,
    'efield_dir': 1,
}


# developers only
neb_sort_tol = 1
sella_eta = 0.002
scale_sigma = 1.0

os.environ['OMP_NUM_THREADS'] = f'{omp}'
profile = AbacusProfile(f'mpirun -np  {mpi} {abacus}')

# reading part
msg = '''
Usage: 
- For using IS and FS: 
    python neb2sella_abacus.py [init_stru] [final_stru] ([format])
- For using existed NEB: 
    python neb2sella_abacus.py [neb_latest.traj]
'''
if len(sys.argv) < 2:
    print(msg)
    sys.exit(1)
elif len(sys.argv) == 2:
    if sys.argv[1] == "-h" or sys.argv[1] == "--help":
        print(msg)
        sys.exit(0)
    else:
        neb_traj = sys.argv[1]
        neb_abacus = read(neb_traj, ":", format="traj")
        atom_init = neb_abacus[0]
        atom_final = neb_abacus[-1]
        assert type(atom_init) == Atoms and type(atom_final) == Atoms, \
        "The input file is not a trajectory file contained Atoms object"
else:
    init_stru = sys.argv[1]
    final_stru = sys.argv[2]
    if len(sys.argv) == 4:
        format = sys.argv[3]
    else:
        format = None # auto detect
    atom_init = read(init_stru, format=format)
    atom_final = read(final_stru, format=format)


# init and final stru
atom_init.calc = Abacus(profile=profile, directory="IS_OPT",
                    **parameters)
init_relax = OPTSolver(atom_init)
init_relax.run(fmax=0.05)
atom_final.calc = Abacus(profile=profile, directory="FS_OPT",
                    **parameters)
final_relax = OPTSolver(atom_final)
final_relax.run(fmax=0.05)

write("init_opted.traj", atom_init, format="traj")
write("final_opted.traj", atom_final, format="traj")

# run neb
print(f'Generating NEB path, number of images: {n_neb_images}, sort_tol: {neb_sort_tol}')
is_pmg = AseAtomsAdaptor.get_structure(atom_init)
fs_pmg = AseAtomsAdaptor.get_structure(atom_final)
path = IDPPSolver.from_endpoints([is_pmg, fs_pmg], n_neb_images, sort_tol=neb_sort_tol)
new_path = path.run(maxiter=5000, tol=1e-5, gtol=1e-3)
# conver path to ase format, and add SinglePointCalculator
ase_path = [i.to_ase_atoms() for i in new_path]
ase_path[0].calc = deepcopy(atom_init.calc)
ase_path[-1].calc = deepcopy(atom_final.calc)
for img in ase_path[1:-1]:
    img.calc = Abacus(profile=profile, directory="NEB",
                    **parameters)

neb = DyNEB(ase_path, 
            climb=climb, dynamic_relaxation=True, fmax=neb_fmax,
            method=neb_algorism, parallel=False, scale_fmax=scale_sigma,
            allow_shared_calculator=True)

traj = Trajectory(neb_log, 'w', neb)
opt = NEBSolver(neb, trajectory=traj)
opt.run(neb_fmax)

# neb displacement to dimer
n_images = NEBTools(ase_path)._guess_nimages()
neb_raw_barrier = max([image.get_potential_energy() for image in ase_path])
fmax = NEBTools(ase_path).get_fmax()
barrier = NEBTools(ase_path).get_barrier()[0]
TS_info = [(ind, image) 
            for ind, image in enumerate(ase_path) 
            if image.get_potential_energy() == neb_raw_barrier][0]
print(f"=== Locate TS in {TS_info[0]} of 0-{n_images-1} images  ===")
print(f"=== NEB Raw Barrier: {neb_raw_barrier:.4f} (eV) ===")
print(f"=== NEB Fmax: {fmax:.4f} (eV/A) ===")
print(f"=== Now Turn to Sella with NEB Information ===")

# sella part
# set cons is optional, ase constraint can be transfered 
ts_guess = TS_info[1].copy()
ts_guess.calc = Abacus(profile=profile, directory="SELLA",
                    **parameters)

dyn = Sella(
    ts_guess,
    trajectory = sella_log,
    eta = sella_eta,
)
dyn.run(fmax=sella_fmax)

# get struc of IS,FS,TS
write("IS_get.cif", atom_init, format="cif")
write("FS_get.cif", atom_final, format="cif")
write("TS_get.cif", ts_guess, format="cif")
write("IS_get.stru", atom_init, format="abacus")
write("FS_get.stru", atom_final, format="abacus")
write("TS_get.stru", ts_guess, format="abacus")

# get energy informations
ene_init = atom_init.get_potential_energy()
ene_final = atom_final.get_potential_energy()
ene_ts = ts_guess.get_potential_energy()
ene_delta = ene_final - ene_init
ene_activa = ene_ts - ene_init
ene_act_rev = ene_ts - ene_final
msg = f'''
==> TS-Search Results <==
- Items      Energy
- IS         {ene_init:.6f}
- FS         {ene_final:.6f}
- TS         {ene_ts:.6f}
- dE         {ene_delta:.6f}
- Ea_f       {ene_activa:.6f}
- Ea_r       {ene_act_rev:.6f}
'''
print(msg)