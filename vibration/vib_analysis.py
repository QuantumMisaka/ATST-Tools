# JamesMisaka in 2023-11-30
# Vibrational analysis from finite displacement by using abacus
# part of ATST-Tools scripts

import os
import numpy as np
from ase.vibrations import Vibrations
from ase.thermochemistry import HarmonicThermo
from ase.io import read, write
from ase.calculators.abacus import Abacus, AbacusProfile
from ase.parallel import world, parprint
from ase.mep.neb import NEBTools
from neb2vib import neb2vib

# reading usage can be changed by user
# usage by neb
neb_traj = read('neb_latest.traj', index=':')
atoms, vib_indices = neb2vib(neb_traj)

# traditional
# stru = "STRU"
# atoms = read(stru)
# vib_indices = [0, 1, 37]

# indices setting for which atoms to be displaced
# vib_indices = [atom.index for atom in atoms if atom.symbol == 'H']

T = 523.15 # K

abacus = "abacus"
mpi = 16
omp = 4
lib_dir = "/lustre/home/2201110432/example/abacus"
pseudo_dir = f"{lib_dir}/PP"
basis_dir = f"{lib_dir}/ORB"
# default pp and basis is supported by ase-abacus interface
pp = {
      'H': 'H_ONCV_PBE-1.0.upf',
      'C': 'C_ONCV_PBE-1.0.upf',
      'O': 'O_ONCV_PBE-1.0.upf',
      'Fe': 'Fe_ONCV_PBE-1.0.upf',
      }
basis = {
         'H': 'H_gga_6au_100Ry_2s1p.orb',
         'C': 'C_gga_7au_100Ry_2s2p1d.orb',
         'O': 'O_gga_7au_100Ry_2s2p1d.orb',
         'Fe': 'Fe_gga_8au_100Ry_4s2p2d1f.orb',
         }
kpts = [3, 1, 2]
parameters = {
    'calculation': 'scf',
    'nspin': 2,
    'xc': 'pbe',
    'ecutwfc': 100,
    'ks_solver': 'genelpa',
    'symmetry': 0,
    'vdw_method': 'none',
    'smearing_method': 'mp',
    'smearing_sigma': 0.002,
    'basis_type': 'lcao',
    'mixing_type': 'broyden',
    'mixing_ndim': 8,
    'scf_thr': 1e-7,
    'scf_nmax': 300,
    'kpts': kpts,
    'pp': pp,
    'basis': basis,
    'pseudo_dir': pseudo_dir,
    'basis_dir': basis_dir,
    'init_chg': 'file',
    'init_wfc': 'atomic',
    'cal_force': 1,
    'cal_stress': 1,
    'out_stru': 1,
    'out_chg': 1,
    'out_mul': 0,
    'out_bandgap': 0,
    'out_wfc_lcao': 0,
    'efield_flag': 1,
    'dip_cor_flag': 1,
    'efield_dir': 1,
}

# developer only
vib_name = 'vib'
delta = 0.01
nfree = 2

def set_calculator(abacus, parameters, mpi=1, omp=1) -> Abacus:
    """Set Abacus calculators"""
    os.environ['OMP_NUM_THREADS'] = f'{omp}'
    profile = AbacusProfile(f"mpirun -np {self.mpi} {self.abacus}")
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
    print("==> Doing Harmonmic Thermodynamic Analysis <==")
    vib_energies = vib.get_energies()
    #real_vib_energies = np.array([energy for energy in vib.get_energies() if energy.imag == 0 and energy.real > 0], dtype=float)
    thermo = HarmonicThermo(vib_energies, ignore_imag_modes=True,)
    entropy = thermo.get_entropy(T)
    free_energy = thermo.get_helmholtz_energy(T)
    print(f"==> Entropy: {entropy:.6e} eV/K <==")
    print(f"==> Free Energy: {free_energy:.6f} eV <==")
    

