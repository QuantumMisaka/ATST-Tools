# JamesMisaka in 2023-11-05
# Run NEB calculation on NEB images by ASE-ABACUS
# part of ASE-NEB-ABACUS scripts

import os
from ase.calculators.abacus import Abacus, AbacusProfile
from ase.optimize import FIRE, BFGS
from ase.mep.neb import NEB, DyNEB # newest ase
from ase.io import read, write
from ase.parallel import world, parprint, paropen
#from pathlib import Path

# setting
mpi = 8
omp = 4
fmax = 0.05  # eV / Ang
neb_optimizer = FIRE # suited for CI-NEB
neb_directory = "NEBrun"
algorism = "improvedtangent" # IT-NEB is recommended
climb = True
dyneb = False  
parallel = True
k = 0.05 # eV/Ang^2, spring constant
init_chain = "init_neb_chain.traj"
abacus = "abacus"
#lib_dir = "/lustre/home/2201110432/example/abacus"
lib_dir = ""
pseudo_dir = f"{lib_dir}"
basis_dir = f"{lib_dir}"
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


class AbacusNEB:
    """Customize Nudged Elastic Band calculation by using ABACUS"""

    def __init__(self, init_chain, parameters, abacus='abacus',
                 dyneb=False, k=0.05, algorism="improvedtangent", 
                 directory='NEB', mpi=1, omp=1, parallel=True, ) -> None:
        """Initialize initial and final states

        init_chain (Atoms object): starting image chain from nem_make.py or other method
        parameters (dict): settings of abacus input parameters
        abacus (str): Abacus executable file. Default: 'abacus'
        dyneb (bool): Using DyNEB method or not, which is efficient in serial NEB calculation. Default: False
        k (float): spring constant for NEB calculation, eV/Ang^2, default: 0.05
        algorism (str): NEB algorism. which can be
        - 'aseneb': standard ase NEB
        - 'improvedtangent' : IT-NEB (recommended by Sobereva)
        - 'eb': climbing image elastic band method (default in AutoNEB)
        - 'spline': 
        - 'string': 
        
        Default: 'improvedtangent'
        
        directory (str): calculator directory name, for parallel calculation {directory}-rank{i} will be the directory name, default: 'NEB'
        mpi (int): number of MPI for abacus calculator
        omp (int): number of OpenMP for abacus calculator
        parallel (bool): parallel calculation setting, default True, if dyneb is True, parallel will be set to False automatically
        """

        self.init_chain = init_chain
        self.dyneb = dyneb
        self.algorism = algorism
        self.abacus = abacus
        self.directory = directory
        self.parameters = parameters
        self.mpi = mpi
        self.omp = omp
        self.parallel = parallel
        self.k = k
        if self.parallel == False:
            parprint("Notice: Parallel calculation is not being used")
            parprint("Set NEB method to DyNEB automatically")
            self.dyneb = True
        elif self.dyneb == True:
            parprint("Notice: Dynamic NEB method is set")
            parprint("Parallel calculation is auto set to False")
            self.parallel = False

    def set_calculator(self):
        """Set Abacus calculators"""
        os.environ['OMP_NUM_THREADS'] = f'{self.omp}'
        profile = AbacusProfile(
            argv=['mpirun', '-np', f'{self.mpi}', self.abacus])
        # add parallel setting 20230926
        if self.parallel:
            out_directory = f"{self.directory}-rank{world.rank}"
        else:
            out_directory = self.directory
        calc = Abacus(profile=profile, directory=out_directory,
                      **self.parameters)
        return calc
    
    
    def set_neb_chain(self, climb=True, fmax=0.05):
        """Set neb_chain, namely:
        1. images defining path from initial to final state
        2. neb object including the images and the implemented NEB method
        
        fmax (float): threshold (unit: eV/Angstrom) of the force convergence
        climb (bool): climbing image NEB method
        interpolate (string): interpolate chain path, 'linear' or 'idpp' or None
        """
        # set calculator
        images = self.init_chain
        for num, image in enumerate(images[1:-1]):
            # only set calc for inter-image
            # we can only use one rank for one inter-image
            if self.parallel:
                if world.rank == num % world.size:
                    image.calc = self.set_calculator()
            # else:
            image.calc = self.set_calculator()
        # setting neb algorism
        if self.algorism in ["aseneb", "improvedtangent", "eb", "spline", "string"]:
            if self.dyneb or not self.parallel :
                # dynamic neb can only be performed by serial
                parprint("----- Running Dynamic NEB -----")
                parprint(f"----- {self.algorism} method is being used -----")
                parprint("----- Default scale_fmax = 1.0 -----")
                # dynamic neb need to read fmax
                neb = DyNEB(images, climb=climb, dynamic_relaxation=True, fmax=fmax
                            method=self.algorism, parallel=False, scale_fmax=1.0, k=self.k)
            else:
                parprint("----- Running ASE-NEB -----")
                parprint(f"----- {self.algorism} method is being used -----")
                if self.parallel:
                    parprint("----- Parallel calculation is being used -----")
                neb = NEB(images, climb=climb, method=self.algorism, 
                          k=self.k, parallel=self.parallel)
        else:
            raise NotImplementedError("NEB algorism not supported! \n Please choose algorism from 'aseneb', 'improvedtangent', 'eb', 'spline', 'string'")
        return neb


    def run(self, optimizer=FIRE, fmax=0.05, climb=True, outfile="neb.traj"):
        """Run Abacus NEB

        optimizer (Optimizer object): defaults to FIRE. BFGS, LBFGS, GPMin, MDMin and QuasiNewton are supported, recommend FIRE method
        fmax (float): threshold (unit: eV/Angstrom) of the force convergence
        climb (bool): climbing image NEB method
        """
        neb = self.set_neb_chain(climb, fmax)
        opt = optimizer(neb, trajectory=outfile)
        opt.run(fmax)
        print("----- NEB calculation finished -----")


if __name__ == "__main__":
# running process
# read initial guessed neb chain
    init_chain = read(init_chain, index=':')
    neb = AbacusNEB(init_chain, parameters=parameters, parallel=parallel,
                    directory=neb_directory, mpi=mpi, omp=omp, abacus=abacus, 
                    algorism=algorism, k=k, dyneb=dyneb)
    neb.run(optimizer=neb_optimizer, climb=climb, fmax=fmax)


