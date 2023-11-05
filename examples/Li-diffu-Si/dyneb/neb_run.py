# JamesMisaka in 2023-10-08
# Run full NEB calculation on NEB images by ASE-ABACUS
# only for Li-diffu-Si example, using dyneb
# part of ASE-NEB-ABACUS scripts

import os 
from ase.calculators.abacus import Abacus, AbacusProfile
from ase.optimize import FIRE, BFGS, QuasiNewton
from ase.mep.neb import NEB, DyNEB # newest ase
from ase.io import read, write
from ase.parallel import world, parprint, paropen
#from pathlib import Path

# setting for NEB
mpi = 4
omp = 2
parallel = False
neb_optimizer = FIRE # suited for CI-NEB
neb_directory = "NEB"
algorism = "improvedtangent" # IT-NEB is recommended
neb_type = "dyneb"
init_chain = "init_neb_chain.traj"
climb = True
fmax = 0.05  # eV / Ang
# setting for calculator
abacus = "abacus"
pseudo_dir = "."
basis_dir = "."
pp = {
    "Li": "Li_ONCV_PBE-1.2.upf",
    "Si": "Si_ONCV_PBE-1.2.upf", }
basis = {
    "Li": "Li_gga_8au_100Ry_4s1p.orb",
    "Si": "Si_gga_8au_100Ry_2s2p1d.orb"}
kpts = [2, 2, 2]
parameters = {
    'calculation': 'scf',
    'xc': 'pbe',
    'ecutwfc': 100,
    'smearing_method': 'gaussian',
    'smearing_sigma': 0.001,
    'basis_type': 'lcao',
    'ks_solver': 'genelpa',
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
    'out_bandgap': 0,
    'efield_flag': 0,
    'dip_cor_flag': 0,
    'efield_dir': 2,
    'efield_pos_max': 0.7
}


class AbacusNEB:
    """Customize Nudged Elastic Band calculation by using ABACUS"""

    def __init__(self, init_chain, parameters, abacus='abacus',
                 neb_type="neb", algorism="improvedtangent", directory='OUT', 
                 mpi=1, omp=1, parallel=False,) -> None:
        """Initialize initial and final states

        init_chain (Atoms object): starting image chain
        parameters (dict): settings of input parameters
        abacus (str): Abacus executable file. Default: 'abacus'
        algorism (str): NEB algorism. which can be
        - 'aseneb': standard ase NEB
        - 'improvedtangent' : IT-NEB (recommended by Sobereva)
        - 'eb': climbing image elastic band method (in AutoNEB)
        - 'spline': 
        - 'string': 
        
        Default: 'improvedtangent'
        
        dyneb (bool): dynamic NEB method. Default: True
        directory (str): work directory
        mpi (int): number of MPI
        omp (int): number of OpenMP
        parallel (bool): parallel calculation setting
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
            self.neb_type = "dyneb"
            

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
                neb = DyNEB(images, climb=climb, fmax=fmax, dynamic_relaxation=True, allow_shared_calculator=True,
                method=self.algorism, parallel=self.parallel, scale_fmax=1.0)
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

## if __name__ == "__main__":
# running process
# read initial guessed neb chain
init_chain = read(init_chain, index=':')

# No atom need to be fixed, pass fix process
# ABACUS fix and magmom information CANNOT be read by ASE-ABACUS now from abacus-out

# do neb calculation parallelly
neb = AbacusNEB(init_chain, parameters=parameters, parallel=parallel,
                directory=neb_directory, mpi=mpi, omp=omp, abacus=abacus, 
                algorism=algorism)
neb.run(optimizer=neb_optimizer, climb=climb, fmax=fmax)


