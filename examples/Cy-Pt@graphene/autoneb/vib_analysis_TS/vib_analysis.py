from ase import Atoms
from ase.optimize import BFGS, FIRE
from ase.vibrations import Vibrations
from ase.io import read, write
from ase.thermochemistry import HarmonicThermo
from ase.calculators.abacus import Abacus, AbacusProfile
from ase.parallel import world, parprint
import os
import numpy as np

#lib_dir = "/home/james/example/"
target = read("neb_latest.traj@-5")
lib_dir = "/lustre/home/2201110432/example/abacus"
mpi = 16
omp = 4

neb_ts = target
vib_indices = list(range(48, 54))
vib_indices += [atom.index for atom in neb_ts if atom.symbol == 'H']
T = 300

pseudo_dir = f"{lib_dir}/PP"
basis_dir = f"{lib_dir}/ORB"
pp = {
      'C':'C_ONCV_PBE-1.0.upf',
      'H':'H_ONCV_PBE-1.0.upf',
      'Pt':'Pt_ONCV_PBE-1.0.upf',
      }
basis = {
         'C': 'C_gga_7au_100Ry_2s2p1d.orb',
         'H': 'H_gga_6au_100Ry_2s1p.orb',
         'Pt': 'Pt_gga_7au_100Ry_4s2p2d1f.orb'
         ,}
kpts = [2, 1, 2]
parameters = {
    'calculation': 'scf',
    'nspin': 2,
    'xc': 'pbe',
    'ecutwfc': 100,
    'ks_solver': 'genelpa',
    'symmetry': 0,
    'vdw_method': 'd3_bj',
    'smearing_method': 'gaussian',
    'smearing_sigma': 0.001,
    'basis_type': 'lcao',
    'mixing_type': 'broyden',
    'scf_thr': 1e-7,
    'scf_nmax': 200,
    'kpts': kpts,
    'pp': pp,
    'basis': basis,
    'pseudo_dir': pseudo_dir,
    'basis_dir': basis_dir,
    'cal_force': 1,
    'cal_stress': 1,
    'out_stru': 1,
    'out_chg': 0,
    'out_mul': 0,
    'out_bandgap': 0,
    'efield_flag': 1,
    'dip_cor_flag': 1,
    'efield_dir': 1,
    'efield_pos_max': 0.0
}

def set_calculator(abacus, parameters, mpi=1, omp=1) -> Abacus:
    """Set Abacus calculators"""
    os.environ['OMP_NUM_THREADS'] = f'{omp}'
    profile = AbacusProfile(
        argv=['mpirun', '-np', f'{mpi}', abacus])
    out_directory = f"OUT-rank{world.rank}"
    calc = Abacus(profile=profile, directory=out_directory,
                **parameters)
    return calc


if __name__ == "__main__":
    neb_ts.calc = set_calculator("abacus", parameters, mpi=mpi, omp=omp)
    vib = Vibrations(neb_ts, indices=vib_indices)
    parprint("==> Running Vibrational Analysis <==")
    vib.run()
    # post-processing
    parprint("==> Get Energy and Frequency <==")
    vib.get_energies()
    vib.get_frequencies()
    parprint("==> Get ZPE <==")
    vib.get_zero_point_energy()
    parprint("==> Writing Mode <==")
    vib.write_mode()
    parprint("==> Summary <==")
    vib.summary()
    # thermochemistry
    print("==> Doing Harmonmic Thermodynamic Analysis <==")
    real_vib_energies = np.array([energy for energy in vib.get_energies() if energy.imag == 0 and energy.real > 0], dtype=float)
    thermo = HarmonicThermo(real_vib_energies)
    entropy = thermo.get_entropy(T)
    free_energy = thermo.get_helmholtz_energy(T)
    print(f"==> Entropy: {entropy:.6e} eV/K <==")
    print(f"==> Free Energy: {free_energy:.6f} eV <==")