from ase.io import read, write
from ase import Atoms
from ase.optimize import BFGS, FIRE, QuasiNewton
from ase.constraints import FixAtoms
from ase.visualize import view
from ase.mep.neb import NEBTools, NEB
from ase.mep.autoneb import AutoNEB
from deepmd_pt.utils.ase_calc import DPCalculator as DP


from ase.mep.autoneb import AutoNEB

model = "FeCHO-dpa2-full.pt"
prefix = "DP-AutoNEB"
n_max = 8
n_simul = 4

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

#atom_init.get_potential_energy()
chain = [atom_init]
for i in range(n_simul):
    image = atom_init.copy()
    image.set_calculator(DP(model=model))
    chain.append(image)
chain.append(atom_final)
neb = NEB(chain)
neb.interpolate(method="idpp")

for num, image in enumerate(chain):
    index = f"{num:03d}"
    write(f"{prefix}{index}.traj", image, format="traj")

def attach_calculators(images: list):
    """Attach calculator to a list of images supplied"""
    for image in images:
        image.calc = DP(model=model)

autoneb = AutoNEB(attach_calculators, prefix=prefix, n_simul=n_simul, n_max=n_max+2, fmax=0.20, climb=True, method="improvedtangent", smooth_curve=True, parallel=False, optimizer=FIRE)
autoneb.run()

