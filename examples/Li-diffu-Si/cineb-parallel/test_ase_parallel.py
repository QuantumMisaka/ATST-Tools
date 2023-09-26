# JamesMisaka in 2023-09-20
# parallel method and AutoNEB needed
import os 
# from ase.neb import NEB, DyNEB, NEBTools # old ase
from ase.mep.neb import NEB, DyNEB, NEBTools # newest ase
# from ase.autoneb import AutoNEB
from ase.mep.autoneb import AutoNEB # newest ase
from ase.calculators.abacus import Abacus, AbacusProfile
from ase.optimize import FIRE, BFGS
from ase.io import read, write
from pathlib import Path
from ase.parallel import world, parprint, paropen

rank = world.rank
size = world.size

for i in range(size):
    if rank == i:
        print(rank, size)