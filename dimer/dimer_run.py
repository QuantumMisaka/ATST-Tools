from ase import Atom, Atoms
from ase.io import Trajectory, read, write
from ase.mep import DimerControl, MinModeAtoms, MinModeTranslate
from ase.optimize import BFGS, FIRE
from ase.calculators.abacus import Abacus, AbacusProfile
import os, sys

import numpy as np

fmax = 0.05
dimer_input_file = 'dimer_init.traj'
init_eigenmode_method = 'displacement'
displacement_input = 'displacement_vector.npy'

# setting for calculator
mpi = 16
omp = 4
abacus = "abacus"
#lib_dir = "/lustre/home/2201110432/example/abacus"
lib_dir = "/home/james/example"
pseudo_dir = f"{lib_dir}/PP"
basis_dir = f"{lib_dir}/ORB"
pp = {
        'H':'H_ONCV_PBE-1.0.upf',
        'Au':'Au_ONCV_PBE-1.0.upf',
      }
basis = {
        'H':'H_gga_6au_100Ry_2s1p.orb',
        'Au':'Au_gga_7au_100Ry_4s2p2d1f.orb',
        }
kpts = [3, 1, 3]
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
    'efield_pos_max': 0.7
}


class AbacusDimer:
    """Customize Dimer calculation by using ABACUS"""
    
    def __init__(self, init_Atoms, parameters, abacus='abacus',
                 mpi=1, omp=1, directory='DIMER', 
                 traj_file='dimer.traj',
                 init_eigenmode_method='displacement',
                 displacement_vector: np.ndarray = None,):
        """Initialize Dimer method by using ASE-ABACUS

        init_Atoms (Atoms object): starting image, can be from every way including NEB result
        parameters (dict): settings of abacus input parameters
        abacus (str): Abacus executable file. Default: 'abacus'
        directory (str): calculator directory name, for parallel calculation {directory}-rank{i} will be the directory name
        mpi (int): number of MPI for abacus calculator
        omp (int): number of OpenMP for abacus calculator
        traj_file (str): trajectory file name for dimer calculation, when running dimer calculation, trajetory will be written to this file, default is 'dimer.traj'
        init_eigenmode_method (str): dimer initial eigenmode method. Choose from 'displacement' and 'gauss'.
        displacement_vector (np.ndarray): displacement vector for dimer initial eigenmode. Only used when init_eigenmode_method is 'displacement'
        """
        self.init_Atoms = init_Atoms
        self.parameters = parameters
        self.abacus = abacus
        self.mpi = mpi
        self.omp = omp
        self.directory = directory
        self.traj_file = traj_file
        self.init_eigenmode_method = init_eigenmode_method
        self.displacement_vector = displacement_vector
        
    def set_calculator(self):
        """Set Abacus calculators"""
        os.environ['OMP_NUM_THREADS'] = f'{self.omp}'
        profile = AbacusProfile(
            argv=['mpirun', '-np', f'{self.mpi}', self.abacus])
        out_directory = self.directory
        calc = Abacus(profile=profile, directory=out_directory,
                        **self.parameters)
        return calc
    
    def set_d_mask_by_displacement(self):
        """set mask by displacement"""
        print("=== Set mask by displacement vector where displacement is [0,0,0] ===")
        d_mask = self.displacement_vector != np.zeros(3)
        d_mask = d_mask[:,0].tolist()
        return d_mask
    
    def set_d_mask_by_constraint(self):
        """set mask by constraint of Atoms"""
        print("=== Set mask by constraint read from init Atoms ===")
        dimer_init = self.init_Atoms
        const_list = dimer_init._get_constraints()[0]
        const_list = const_list.todict()['kwargs']['indices']
        d_mask = [True] * len(dimer_init)
        for ind in const_list:
            d_mask[ind] = False
        return d_mask
        
    def run(self, fmax=0.05):
        """run dimer calculation workflow"""
        dimer_init = self.init_Atoms
        dimer_init.calc = self.set_calculator()
        dimer_traj = Trajectory(self.traj_file, 'w', dimer_init)
        if self.init_eigenmode_method == "displacement":
            # d_mask = self.set_d_mask_by_displacement() # not the best
            d_mask = self.set_d_mask_by_constraint()
            d_control = DimerControl(initial_eigenmode_method=self.init_eigenmode_method, 
                                    displacement_method="vector", 
                                    mask=d_mask)
            d_atoms = MinModeAtoms(dimer_init, d_control)
            d_atoms.displace(displacement_vector=self.displacement_vector)
        elif self.init_eigenmode_method == "gauss":
            # leave a way for random displacement
            d_mask = self.set_d_mask_by_constraint()
            d_control = DimerControl(initial_eigenmode_method=self.init_eigenmode_method, 
                                    mask=d_mask)
            d_atoms = MinModeAtoms(dimer_init, d_control)
        else:
            raise ValueError("init_eigenmode_method must be displacement or gauss")
        dimer_relax = MinModeTranslate(d_atoms, trajectory=dimer_traj)
        dimer_relax.run(fmax=fmax)
        
if __name__ == "__main__":
# running process
# read initial guessed neb chain
    dimer_init = read(dimer_input_file)
    displacement_vector = np.load(displacement_input)
    dimer = AbacusDimer(dimer_init, parameters, abacus=abacus,
                        mpi=mpi, omp=omp, 
                        init_eigenmode_method=init_eigenmode_method,
                        displacement_vector=displacement_vector)
    dimer.run(fmax=fmax)
