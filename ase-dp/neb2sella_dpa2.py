import numpy as np
import os

from ase.io import read, write, Trajectory
from ase import Atoms
from ase.optimize import BFGS, FIRE, QuasiNewton
from ase.constraints import FixAtoms
from ase.visualize import view
from ase.mep.neb import NEBTools, NEB, DyNEB
from ase.mep.autoneb import AutoNEB
from sella import Sella, Constraints

from deepmd_pt.utils.ase_calc import DPCalculator as DP


from ase.mep.autoneb import AutoNEB

model = "FeCHO-dpa2-full.pt"
prefix = "DP-AutoNEB"
n_max = 8
neb_fmax = 1.00  # neb should be rough
sella_fmax = 0.05 # dimer use neb guess
climb = True
scale_fmax = 1.0 # use dyneb to reduce message far from TS
omp = 16
neb_algorism = "improvedtangent"
neb_traj = "neb_images.traj"

os.environ['OMP_NUM_THREADS'] = "omp"

# init and final stru
atom_init = read("./C-C-Fe5C2_510/data/IS/STRU")
# atom_init.calc = DP(model=model)
# init_relax = BFGS(atom_init)
# init_relax.run(fmax=0.05)
print(atom_init.get_potential_energy())
atom_final = read("./C-C-Fe5C2_510/data/FS/STRU")
atom_final.calc = DP(model=model)
final_relax = BFGS(atom_final)
final_relax.run(fmax=0.05)

write("init_opted.traj", atom_init, format="traj")
write("final_opted.traj", atom_final, format="traj")

# run neb
images = [atom_init]
for i in range(n_max):
    image = atom_init.copy()
    image.set_calculator(DP(model=model))
    images.append(image)
images.append(atom_final)
neb = DyNEB(images, 
            climb=climb, dynamic_relaxation=True, fmax=neb_fmax,
            method=neb_algorism, parallel=False, scale_fmax=scale_fmax)
neb.interpolate(method="idpp")

traj = Trajectory(neb_traj, 'w', neb)
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
print(f"=== Now Turn to Dimer with NEB Information ===")



# sella part

# set cons is optional
cons = Constraints()
cons.fix_translation(ts_neb._get_constraints()[0].get_indices())

dimer_init = TS_info[1].copy()
init_eigenmode_method = "displacement"
dimer = DPDimer(dimer_init, model=model,
                        omp=omp, 
                        init_eigenmode_method=init_eigenmode_method,
                        displacement_vector=displacement_vector)
dimer.run(fmax=sella_fmax)