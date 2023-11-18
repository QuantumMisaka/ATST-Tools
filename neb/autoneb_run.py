# JamesMisaka in 2023-11-05
# Run AutoNEB calculation on NEB images by ASE-ABACUS
# part of ASE-NEB-ABACUS scripts

import os
from ase.calculators.abacus import Abacus, AbacusProfile
from ase.optimize import FIRE, BFGS
from ase.mep.autoneb import AutoNEB
from ase.io import read, write
from ase.parallel import world, parprint, paropen
#from pathlib import Path

# setting
mpi = 16
omp = 4
neb_optimizer = FIRE # suited for CI-NEB
neb_directory = "AutoNEBrun"
# algorism = 'eb' # default for AutoNEB
algorism = "improvedtangent" # IT-NEB
init_chain = "init_neb_chain.traj"
climb = True
fmax = [0.20, 0.05]  # eV / Ang, 2 fmax for others and last CI-NEB in all-images
n_simul = world.size # only for autoneb, number of simultaneous calculation
n_images = 10 # only for autoneb, max number of all image, which should be reached
smooth_curve = False # True to do more neb step to smooth the curve.
k = 0.05 # eV/Ang^2, force constant of spring, 0.05 is from VTST-Tools
abacus = "abacus"
#lib_dir = "/lustre/home/2201110432/example/abacus"
lib_dir = ""
pseudo_dir = f"{lib_dir}/"
basis_dir = f"{lib_dir}/"
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


class AbacusAutoNEB:
    """Customize AutoNEB calculation by using ABACUS"""

    def __init__(self, init_chain, parameters, 
                 abacus='abacus',  prefix="run_autoneb", 
                 n_simul=1, n_max=10, k=0.05,
                 algorism="improvedtangent", 
                 directory='AutoNEBrun', mpi=1, omp=1, parallel=True, ):
        """Initialize initial and final states

        init_chain (Atoms object): starting image chain from nem_make.py or other method, can only include Initial and Final states, autoneb will generate inter-images automactically. 
        parameters (dict): settings of abacus input parameters
        abacus (str): Abacus executable file. Default: 'abacus'
        algorism (str): NEB algorism. which can be
        - 'aseneb': standard ase NEB
        - 'improvedtangent' : IT-NEB (recommended by Sobereva)
        - 'eb': climbing image elastic band method (default in AutoNEB)
        
        Default: 'improvedtangent'
        
        prefix (str): prefix for AutoNEB output files, default 'autoneb'
        n_simul (int): number of simultaneous calculation, default use world.size, read from ase.parallel
        n_max (int): max number of all image in NEB band, default 10. The max number will not be reached if the convergence is reached before that.
        directory (str): calculator directory name, for parallel calculation {directory}-rank{i} will be the directory name
        mpi (int): number of MPI for abacus calculator
        omp (int): number of OpenMP for abacus calculator
        parallel (bool): parallel calculation setting, default True
        """

        self.init_chain = init_chain
        self.algorism = algorism
        self.abacus = abacus
        self.directory = directory
        self.parameters = parameters
        self.prefix = prefix
        self.mpi = mpi
        self.omp = omp
        self.k = k
        self.parallel = parallel
        parprint("Notice: AutoNEB method is set")
        if (n_simul > 0) and (n_max >= n_simul):
            parprint(f"You manually set n_simul = {n_simul}, n_max = {n_max}", )
            self.n_simul = n_simul
            self.n_max = n_max
        else:
            raise ValueError("You must set n_simul > 0 and n_max >= n_simul numbers for AutoNEB")
        
            
    def set_calculator(self):
        """Set Abacus calculators"""
        os.environ['OMP_NUM_THREADS'] = f'{self.omp}'
        profile = AbacusProfile(
            argv=['mpirun', '-np', f'{self.mpi}', self.abacus])
        if self.parallel:
            out_directory = f"{self.directory}-rank{world.rank}"
        else:
            out_directory = self.directory
        calc = Abacus(profile=profile, directory=out_directory,
                    **self.parameters)
        return calc
    
    
    def attach_calculators(self, images: list):
        """Attach calculator to a list of images supplied"""
        for num, image in enumerate(images):
            image.calc = self.set_calculator()
        
    
    def run(self, optimizer=FIRE, fmax=0.05, climb=True, smooth_curve=False):
        """Run Abacus AutoNEB

        optimizer (Optimizer object): defaults to FIRE. BFGS and FIRE is only used
        fmax (float): threshold (unit: eV/Angstrom) of the force convergence
        climb (bool): climbing image NEB method
        """
        parprint("----- Running AutoNEB -----")
        parprint(f"----- {self.algorism} method is being used -----")
        for num, image in enumerate(self.init_chain):
            index = f"{num:03d}"
            write(f"{self.prefix}{index}.traj", image, format="traj")
        autoneb = AutoNEB(self.attach_calculators, self.prefix, self.n_simul, 
                          self.n_max, fmax=fmax, climb=climb, k = self.k,
                            method=self.algorism, parallel=self.parallel, 
                            optimizer=optimizer, smooth_curve=smooth_curve)
        autoneb.run()
        parprint("----- AutoNEB calculation finished -----")


if __name__ == "__main__": 
# running process
# read initial guessed neb chain
    init_chain = read(init_chain, index=':')
    neb = AbacusAutoNEB(init_chain, parameters, algorism=algorism, 
                        directory=neb_directory, k=k,
                        n_simul=n_simul, n_max=n_images, 
                        abacus=abacus,  mpi=mpi, omp=omp, )
    neb.run(optimizer=neb_optimizer, climb=climb, 
                fmax=fmax, smooth_curve=smooth_curve)