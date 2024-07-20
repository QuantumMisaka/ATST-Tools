from ase.io import read, write
from ase import Atoms
from ase.optimize import BFGS, FIRE, QuasiNewton
from ase.io import Trajectory, read, write
import os
import numpy as np

from sella import Sella, Constraints

from deepmd_pt.utils.ase_calc import DPCalculator as DP

input_file = 'STRU'
model = "FeCHO-dpa2-full.pt"
fmax = 0.05 # sella use neb guess
sella_log = "sella_images.traj"
# setting for calculator
omp = 16



if __name__ == "__main__":
    # running process
    stru = read(input_file)
    stru.calc = DP(model=model)
    dyn = Sella(
    stru,
    trajectory=sella_log,
    )
    dyn.run(fmax=fmax)
    # output ts stru
    write("sella_opted.cif", stru, format="cif")
    write("sella_opted.stru", stru, format="abacus")
    