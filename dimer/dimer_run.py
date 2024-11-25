# JamesMisaka in 2024-09-16
# Run Dimer calculation based on initial guess by ASE-ABACUS
# part of ATST-Tools scripts

from ase.io import Trajectory, read, write
import numpy as np
from abacus_dimer import AbacusDimer

fmax = 0.05
dimer_input_file = 'dimer_init.traj' # "STRU"
init_eigenmode_method = 'displacement'
displacement_input = 'displacement_vector.npy'

# setting for calculator
# abacus calculator setting
abacus = "abacus"
mpi=16
omp=4
moving_atoms_ind = None
lib_dir = "/lustre/home/2201110432/example/abacus"
#lib_dir = "/data/home/liuzq/example"
pseudo_dir = f"{lib_dir}/PP"
basis_dir = f"{lib_dir}/ORB"
pp = {
        "H":"H_ONCV_PBE-1.0.upf",
        "C":"C_ONCV_PBE-1.0.upf",
        "O":"O_ONCV_PBE-1.0.upf",
        "Fe":"Fe_ONCV_PBE-1.0.upf",
      }
basis = {
        "H":"H_gga_6au_100Ry_2s1p.orb",
        "C":"C_gga_7au_100Ry_2s2p1d.orb",
        "O":"O_gga_7au_100Ry_2s2p1d.orb",
        "Fe":"Fe_gga_8au_100Ry_4s2p2d1f.orb",
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
    'mixing_ndim': 20,
    'scf_thr': 1e-7,
    'scf_nmax': 100,
    'kpts': kpts,
    'pp': pp,
    'basis': basis,
    'pseudo_dir': pseudo_dir,
    'basis_dir': basis_dir,
    'cal_force': 1,
    'cal_stress': 1,
    'init_wfc': 'atomic',
    'init_chg': 'auto',
    'out_stru': 1,
    'out_chg': 0,
    'out_mul': 1,
    'out_wfc_lcao': 0,
    'out_bandgap': 0,
    'efield_flag': 1,
    'dip_cor_flag': 1,
    'efield_dir': 1,
}

        
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

    # output TS stru file
    write("TS_dimer.stru", dimer_init, format="abacus")
    write("TS_dimer.cif", dimer_init, format="cif")