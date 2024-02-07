from ase.io import read, write, Trajectory
from ase import Atoms
from ase.optimize import BFGS, FIRE, QuasiNewton
from ase.constraints import FixAtoms
from ase.visualize import view
from ase.mep.neb import NEBTools, NEB, DyNEB
from ase.mep.autoneb import AutoNEB
from deepmd_pt.utils.ase_calc import DPCalculator as DP


from ase.mep.autoneb import AutoNEB

model = "FeCHO-dpa2-full.pt"
prefix = "DP-AutoNEB"
n_max = 8
n_simul = 4

neb_fmax = 0.20
dimer_fmax = 0.05
climb = True
scale_fmax = 1.0
neb_algorism = "improvedtangent"
neb_traj = "neb_images.traj"

# init and final stru
atom_init = read("./C-C-Fe5C2_510/data/IS/STRU")
atom_final = read("./C-C-Fe5C2_510/data/FS/STRU")
atom_init.calc = DP(model=model)
atom_final.calc = DP(model=model)
init_relax = BFGS(atom_init)
final_relax = BFGS(atom_final)
init_relax.run(fmax=0.05)
final_relax.run(fmax=0.05)

write("init_opted.traj", atom_init, format="traj")
write("final_opted.traj", atom_final, format="traj")

# run neb
images = [atom_init]
for i in range(n_simul):
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
n_images = NEBTools(neb_traj)._guess_nimages()


