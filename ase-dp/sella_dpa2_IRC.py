from ase.io import read, write
from ase import Atoms
from ase.optimize import BFGS, FIRE, QuasiNewton
from ase.io import Trajectory, read, write
import os
import numpy as np

from sella import Sella, Constraints, IRC

from deepmd_pt.utils.ase_calc import DPCalculator as DP

input_file = 'STRU'
model = "FeCHO-dpa2-full.pt"
fmax = 0.05 # sella use neb guess
irc_log = "sella_IRC.traj"
steps = 1000 # IRC steps
dx = 0.1, # IRC step size
# setting for calculator
omp = 16



if __name__ == "__main__":
    # running process
    stru = read(input_file)
    stru.calc = DP(model=model)
    irc = IRC(stru, irc_log, dx=dx,)
    irc.run(fmax, steps, direction='forward')
    irc.run(fmax, steps, direction='backward')
    
