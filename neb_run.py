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

# setting for NEB
mpi = 16
omp = 2
neb_optimizer = FIRE # suited for CI-NEB
neb_directory = "NEBrun"
algorism = "improvedtangent" # IT-NEB is recommended
neb_type = "neb"  # neb, dyneb, but not autoneb
init_chain = "init_neb_chain.traj"
climb = True
fmax = 0.05  # eV / Ang
parallel = True

# setting for calculator
abacus = "abacus"
lib_dir = "/lustre/home/2201110432/example/abacus"
#lib_dir = ""
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
    'symmetry': 1,
    'symmetry_autoclose': 1,
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
                 neb_type="neb", algorism="improvedtangent", directory='NEB', 
                 mpi=1, omp=1, parallel=True, n_simul_auto=0, n_max_auto=0) -> None:
        """Initialize initial and final states

        init_chain (Atoms object): starting image chain from nem_make.py or other method
        parameters (dict): settings of abacus input parameters
        abacus (str): Abacus executable file. Default: 'abacus'
        algorism (str): NEB algorism. which can be
        - 'aseneb': standard ase NEB
        - 'improvedtangent' : IT-NEB (recommended by Sobereva)
        - 'eb': climbing image elastic band method (default in AutoNEB)
        - 'spline': 
        - 'string': 
        
        Default: 'improvedtangent'
        
        neb_type (str): NEB method.  Choose from 'neb', 'dyneb'. and 'autoneb' should use AbacusAutoNEB class
        directory (str): calculator directory name, for parallel calculation {directory}-rank{i} will be the directory name
        mpi (int): number of MPI for abacus calculator
        omp (int): number of OpenMP for abacus calculator
        parallel (bool): parallel calculation setting, default True, if neb_type set to 'dyneb', parallel will be set to False automatically
        """

        self.init_chain = init_chain
        self.neb_type = neb_type
        self.algorism = algorism
        self.abacus = abacus
        self.directory = directory
        self.parameters = parameters
        self.mpi = mpi
        self.omp = omp
        self.parallel = parallel
        if self.parallel == False:
            parprint("Notice: Parallel calculation is not being used")
            parprint("Set NEB method to DyNEB automatically")
            self.neb_type = 'dyneb'
        elif self.algorism == 'dyneb':
            parprint("Notice: Dynamic NEB method is set")
            parprint("Parallel calculation is auto set to False")
            self.parallel = False
        if self.neb_type == 'autoneb':
            parprint("Notice: AutoNEB method is set")
            if (n_simul_auto > 0) and (n_max_auto >= n_simul_auto):
                parprint(f"You manually set n_simul = {n_simul_auto}, n_max = {n_max_auto}", )
                self.n_simul_auto = n_simul_auto
                self.n_max_auto = n_max_auto
            else:
                raise ValueError("You must set n_simul > 0 and n_max >= n_simul numbers for AutoNEB")

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
    
    
    def set_neb_chain(self, fmax=0.05, climb=True,):
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
            if self.neb_type == "dyneb" or not self.parallel :
                # dynamic neb can only be performed by serial
                parprint("----- Running Dynamic NEB -----")
                parprint(f"----- {self.algorism} method is being used -----")
                parprint("----- Default scale_fmax = 1.0 -----")
                neb = DyNEB(images, climb=climb, fmax=fmax, dynamic_relaxation=True, 
                            method=self.algorism, parallel=False, scale_fmax=1.0)
            elif self.neb_type == "neb":
                parprint("----- Running ASE-NEB -----")
                parprint(f"----- {self.algorism} method is being used -----")
                if self.parallel:
                    parprint("----- Parallel calculation is being used -----")
                neb = NEB(images, climb=climb, method=self.algorism, parallel=self.parallel)
            else:
                raise NotImplementedError("----- neb_type NOT SUPPORTED ! -----")
        else:
            raise NotImplementedError("NEB algorism not supported! \n Please choose algorism from 'aseneb', 'improvedtangent', 'eb', 'spline', 'string'")
        return neb


    def run(self, optimizer=FIRE, fmax=0.05, climb=True, outfile="neb.traj"):
        """Run Abacus NEB

        optimizer (Optimizer object): defaults to FIRE. BFGS, LBFGS, GPMin, MDMin and QuasiNewton are supported, recommend FIRE method
        fmax (float): threshold (unit: eV/Angstrom) of the force convergence
        climb (bool): climbing image NEB method
        """
        neb = self.set_neb_chain(fmax, climb)
        opt = optimizer(neb, trajectory=outfile)
        opt.run(fmax)
        print("----- NEB calculation finished -----")


if __name__ == "__main__":
# running process
# read initial guessed neb chain
    init_chain = read(init_chain, index=':')
    neb = AbacusNEB(init_chain, parameters=parameters, parallel=parallel,
                    directory=neb_directory, mpi=mpi, omp=omp, abacus=abacus, 
                    algorism=algorism)
    neb.run(optimizer=neb_optimizer, climb=climb, fmax=fmax)


