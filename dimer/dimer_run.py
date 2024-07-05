# JamesMisaka in 2023-11-27
# Run Dimer calculation based on initial guess by ASE-ABACUS
# part of ATST-Tools scripts

from ase.io import Trajectory, read, write
import os, sys
import numpy as np
from abacus_dimer import AbacusDimer

fmax = 0.05
dimer_input_file = 'dimer_init.traj'
init_eigenmode_method = 'displacement'
displacement_input = 'displacement_vector.npy'

# setting for calculator
mpi = 16
omp = 4
abacus = "abacus"
moving_atoms_ind = None
#lib_dir = "/lustre/home/2201110432/example/abacus"
lib_dir = ""
pseudo_dir = f"{lib_dir}/PP"
basis_dir = f"{lib_dir}/ORB"
# default pp and basis is supported by ase-abacus interface
# pp = {
#       'C':'C_ONCV_PBE-1.0.upf',
#       'H':'H_ONCV_PBE-1.0.upf',
#       'Pt':'Pt_ONCV_PBE-1.0.upf',
#       }
# basis = {
#          'C': 'C_gga_7au_100Ry_2s2p1d.orb',
#          'H': 'H_gga_6au_100Ry_2s1p.orb',
#          'Pt': 'Pt_gga_7au_100Ry_4s2p2d1f.orb'
#          ,}
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
    'scf_thr': 1e-6,
    'scf_nmax': 100,
    'kpts': kpts,
    'pseudo_dir': pseudo_dir,
    'basis_dir': basis_dir,
    'init_chg': 'atomic',
    'init_wfc': 'atomic',
    'cal_force': 1,
    'cal_stress': 1,
    'out_stru': 1,
    'out_chg': 0,
    'out_mul': 0,
    'out_bandgap': 0,
    'out_wfc_lcao': 0,
    'efield_flag': 1,
    'dip_cor_flag': 1,
    'efield_dir': 1,
}
    # 'pp': pp,
    # 'basis': basis,
        
if __name__ == "__main__":
# running process
# read initial guessed neb chain
    dimer_init = read(dimer_input_file)
    displacement_vector = np.load(displacement_input)
    dimer = AbacusDimer(dimer_init, parameters, abacus=abacus,
                        mpi=mpi, omp=omp, 
                        init_eigenmode_method=init_eigenmode_method,
                        displacement_vector=displacement_vector)
    dimer.run(fmax=fmax, moving_atoms_ind=moving_atoms_ind)
