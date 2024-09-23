# JamesMisaka in 2024-01-16
# Thermodynamic analysis for molecular (in a box)
# part of ATST-Tools scripts

import os
import numpy as np
from ase.vibrations import Vibrations
from ase.thermochemistry import IdealGasThermo
from ase.io import read, write
from ase.calculators.abacus import Abacus, AbacusProfile
from ase.parallel import world, parprint
from ase import units
from deepmd.calculator import DP

# setting
model = "FeCHO-dpa2-full.pt"
stru = "STRU"
atoms = read(stru)
# indices setting for which atoms to be displaced
#vib_indices = [atom.index for atom in atoms if atom.symbol == 'H']
vib_indices = None
T = 523.15 # K
P = 1E5 # Pa
geometry = 'linear'
symmetrynumber = 1

vib_name = 'vib'
delta = 0.01
nfree = 2

if __name__ == "__main__":
    print("==> Starting Vibrational Analysis <==")

    atoms.calc = DP(model=model)
    vib = Vibrations(atoms, indices=vib_indices, 
                    name=vib_name, delta=delta, nfree=nfree)

    print("==> Running Vibrational Analysis <==")
    vib.run()
    # post-processing
    print("==> Done !!! <==")
    print(f"==> All force cache will be in {vib_name} directory <==")
    print("==> Vibrational Analysis Summary <==")
    vib.summary()
    print("==> Writing All Mode Trajectory <==")
    vib.write_mode()
    # thermochemistry
    print("==> Doing Ideal-Gas Thermodynamic Analysis <==")
    # initias = (atoms.get_moments_of_inertia() * units._amu /
    #                     (10.0**10)**2)
    # print(f"==> Initial Moments of Inertia for {geometry} {atoms.symbols} molecular: {initias} kg*m^2 <==")
    gasthermo = IdealGasThermo(vib.get_energies(),
                           geometry=geometry,
                           atoms=atoms,
                           ignore_imag_modes=True,
                           symmetrynumber=symmetrynumber,
                           spin=0,)
    entropy = gasthermo.get_entropy(T,P)
    free_energy = gasthermo.get_gibbs_energy(T,P)
    print(f"==> Entropy: {entropy:.6e} eV/K <==")
    print(f"==> Gibbs Free Energy: {free_energy:.6f} eV <==")
    
    

