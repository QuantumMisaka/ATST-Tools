# read STRU file, run optimization jobs
import os
import sys
from ase.optimize import QuasiNewton, BFGS, LBFGS, FIRE, GPMin, MDMin
from ase.io import read, write
from ase.io import Trajectory
from deepmd_pt.utils.ase_calc import DPCalculator as DP

# setting
model = "FeCHO-dpa2-full.pt"
optimizer = BFGS
omp = 4
out_traj = "relax.traj"

if len(sys.argv) > 1:
    stru = read(sys.argv[1], format='abacus')
else:
    stru = read('STRU', format='abacus')

# running
os.environ['OMP_NUM_THREADS'] = f'{omp}'
stru.calc = DP(model=model)
qn = optimizer(stru, trajectory=out_traj)
qn.run(fmax=0.05)
# print-out
print(f"Final energy: {stru.get_potential_energy()} eV")
write("STRU_latest", stru, format='abacus')
write("STRU_latest.cif", stru, format='cif')
