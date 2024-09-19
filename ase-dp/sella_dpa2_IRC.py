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
dx = 0.1 # IRC step size
# setting for calculator
omp = 16
# for developers
sella_eta = 2e-4

if __name__ == "__main__":
    # running process
    irc_traj = Trajectory(irc_log, 'w')
    stru = read(input_file)
    stru.calc = DP(model=model)
    irc = IRC(stru, trajectory=irc_traj, dx=dx, eta=sella_eta)
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
