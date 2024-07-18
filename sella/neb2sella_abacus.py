# need test
import numpy as np
import os, sys

from ase.io import read, write, Trajectory
from ase import Atoms
from ase.optimize import BFGS, FIRE, QuasiNewton
from ase.constraints import FixAtoms
from ase.visualize import view
from ase.mep.neb import NEBTools, NEB, DyNEB
from ase.mep.autoneb import AutoNEB
from sella import Sella, Constraints
from ase.calculators.abacus import Abacus, AbacusProfile

# from deepmd_pt.utils.ase_calc import DPCalculator as DP


n_max = 8
neb_fmax = 0.80  # neb should be rough
sella_fmax = 0.05 # dimer use neb guess
climb = True
scale_fmax = 1.0 # use dyneb to reduce message far from TS
neb_algorism = "improvedtangent"
neb_log = "neb_images.traj"
sella_log = "sella_images.traj"
abacus = "abacus"

# abacus setting
mpi = 16
omp = 4
abacus = "abacus"
lib_dir = "/lustre/home/2201110432/example/abacus"
#lib_dir = "/data/home/liuzq/example"
pseudo_dir = f"{lib_dir}/PP"
basis_dir = f"{lib_dir}/ORB"
properties = ["energy", "forces", "stress"]
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
    'ks_solver': 'genelpa',
    'vdw_method': 'none',
    'nspin' : 2,
    'xc': 'pbe',
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
    'mixing_ndim': 8,
    'scf_thr': 1e-6,
    'kpts': kpts,
    'init_wfc': 'atomic',
    'init_chg': 'atomic',
    'cal_force': 1,
    'cal_stress': 1,
    'out_stru': 1,
    'out_chg': 0,
    'out_mul': 1,
    'out_bandgap': 1,
    'out_wfc_lcao': 0,
    'efield_flag': 1,
    'dip_cor_flag': 1,
    'efield_dir': 1,
}

os.environ['OMP_NUM_THREADS'] = f'{omp}'
profile = AbacusProfile(f'mpirun -np  {mpi} {abacus}')


# reading part
msg = '''
Usage: 
- For using IS and FS: 
    python neb2dimer_dpa2.py [init_stru] [final_stru] ([format])
- For using existing NEB: 
    python neb2dimer_dpa2.py [neb_latest.traj]
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
atom_init.calc = Abacus(profile=profile, directory="IS",
                    **parameters)
# init_relax = BFGS(atom_init)
# init_relax.run(fmax=0.05)
print(atom_init.get_potential_energy())
atom_final.calc = Abacus(profile=profile, directory="FS",
                    **parameters)
final_relax = BFGS(atom_final)
final_relax.run(fmax=0.05)

write("init_opted.traj", atom_init, format="traj")
write("final_opted.traj", atom_final, format="traj")

# run neb
images = [atom_init]
for i in range(n_max):
    image = atom_init.copy()
    image.set_calculator(Abacus(profile=profile, directory="NEB",
                    **parameters))
    images.append(image)
images.append(atom_final)
neb = DyNEB(images, 
            climb=climb, dynamic_relaxation=True, fmax=neb_fmax,
            method=neb_algorism, parallel=False, scale_fmax=scale_fmax)
neb.interpolate(method="idpp")

traj = Trajectory(neb_log, 'w', neb)
opt = FIRE(neb, trajectory=traj)
opt.run(neb_fmax)

# neb displacement to dimer
n_images = NEBTools(images)._guess_nimages()
neb_raw_barrier = max([image.get_potential_energy() for image in images])
fmax = NEBTools(images).get_fmax()
barrier = NEBTools(images).get_barrier()[0]
TS_info = [(ind, image) 
            for ind, image in enumerate(images) 
            if image.get_potential_energy() == neb_raw_barrier][0]
print(f"=== Locate TS in {TS_info[0]} of 0-{n_images-1} images  ===")
print(f"=== NEB Raw Barrier: {neb_raw_barrier:.4f} (eV) ===")
print(f"=== NEB Fmax: {fmax:.4f} (eV/A) ===")
print(f"=== Now Turn to Sella with NEB Information ===")

# para for neb2dimer
# step_after_TS = 1
# norm_vector = 0.01
#out_vec = 'displacement_vector.npy',

# ind_before_TS = TS_info[0] - step_before_TS
# ind_after_TS = TS_info[0] + step_after_TS
# img_before = images[ind_before_TS]
# img_after = images[ind_after_TS]
# image_vector = (img_before.positions - img_after.positions)
# modulo_norm = np.linalg.norm(image_vector) / norm_vector
# displacement_vector = image_vector / modulo_norm
# print(f"=== Displacement vector generated by {ind_before_TS} and {ind_after_TS} images of NEB chain ===")
# print(f"=== Which is normalized to {norm_vector} length ! ===")

# sella part
# set cons is optional
# there are some problems during setting constraints
ts_guess = TS_info[1].copy()
ts_guess.calc = Abacus(profile=profile, directory="SELLA",
                    **parameters)


dyn = Sella(
    ts_guess,
    trajectory=sella_log,
)
dyn.run(fmax=sella_fmax)