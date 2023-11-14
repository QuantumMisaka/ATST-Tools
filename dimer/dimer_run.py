from ase import Atom, Atoms
from ase.io import Trajectory, read, write
from ase.mep import DimerControl, MinModeAtoms, MinModeTranslate
from ase.optimize import BFGS, FIRE

import numpy as np

dimer_input_file = 'dimer_init.traj'
fmax = 0.05
init_eigenmode_method = 'displacement'
displacement_vector = np.zeros(3)



