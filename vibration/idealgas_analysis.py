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

abacus = "abacus"
mpi = 4
omp = 4
lib_dir = "/home/james/example/abacus"
pseudo_dir = f"{lib_dir}/PP"
basis_dir = f"{lib_dir}/ORB"
# default pp and basis is supported by ase-abacus interface
kpts = [1, 1, 1]
parameters = {
    'calculation': 'scf',
    'nspin': 2,
    'xc': 'pbe',
    'ecutwfc': 100,
    'ks_solver': 'genelpa',
    'symmetry': 0,
    'vdw_method': 'none',
    'gamma_only': 1,
    'smearing_method': 'gau',
    'smearing_sigma': 0.001,
    'basis_type': 'lcao',
    'mixing_type': 'broyden',
    'scf_thr': 1e-7,
    'scf_nmax': 300,
    'kpts': kpts,
    'pseudo_dir': pseudo_dir,
    'basis_dir': basis_dir,
    'init_chg': 'atomic',
    'init_wfc': 'file',
    'cal_force': 1,
    'cal_stress': 1,
    'out_stru': 1,
    'out_chg': 1,
    'out_mul': 0,
    'out_bandgap': 0,
    'out_wfc_lcao': 1,
}


def set_calculator(abacus, parameters, mpi=1, omp=1) -> Abacus:
    """Set Abacus calculators"""
    os.environ['OMP_NUM_THREADS'] = f'{omp}'
    profile = AbacusProfile(
        argv=['mpirun', '-np', f'{mpi}', abacus])
    out_directory = f"SCF-rank{world.rank}"
    calc = Abacus(profile=profile, directory=out_directory,
                **parameters)
    return calc

if __name__ == "__main__":
    print("==> Starting Vibrational Analysis <==")

    atoms.calc = set_calculator(abacus, parameters, mpi=mpi, omp=omp)

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
    
    

