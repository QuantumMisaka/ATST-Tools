# JamesMisaka in 2023-09-20
# parallel method and AutoNEB needed
import os 
# from ase.neb import NEB, DyNEB, NEBTools # old ase
from ase.mep.neb import NEB, DyNEB, NEBTools # newest ase
# from ase.autoneb import AutoNEB
from ase.mep.autoneb import AutoNEB # newest ase
from ase.calculators.abacus import Abacus, AbacusProfile
from ase.optimize import FIRE, BFGS
from ase.io import read, write
from pathlib import Path
from ase.parallel import world

# os.environ['ABACUS_PP_PATH'] = '/home/james/example/PP'
# os.environ['ABACUS_ORBITAL_PATH'] = '/home/james/example/ORB'
# print(os.environ.get('ABACUS_PP_PATH'))
# print(os.environ.get('ABACUS_ORBITAL_PATH'))


class AbacusNEB:
    """Customize Abacus workflow for Nudged Elastic Band calculation"""

    def __init__(self, initial, final, parameters, abacus='abacus',
                 algorism="improvedtangent", dyneb=True, directory='OUT', 
                 mpi=1, omp=1, guess=[], n_max=8, parallel=False,) -> None:
        """Initialize initial and final states

        initial (Atoms object): relaxed starting image
        final (Atoms object): relaxed ending image
        n_max (int): the number of images to interpolate between the initial and final images.
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
        guess (List(Atoms)): intermediate guesses. If set, not need to set n_max
        parallel (bool): parallel calculation setting, will be false if use DyNEB
        """

        self.initial = initial
        self.final = final
        self.n_max = n_max
        self.parameters = parameters
        self.algorism = algorism
        self.dyneb = dyneb
        self.abacus = abacus
        self.directory = directory
        self.mpi = mpi
        self.omp = omp
        self.guess = guess
        if self.dyneb == True:
            print("Notice: Dynamic NEB should be efficient in serial calculation")
        self.parallel = parallel
            

    def set_calculator(self):
        """Set Abacus calculators"""
        os.environ['OMP_NUM_THREADS'] = f'{self.omp}'
        profile = AbacusProfile(
            argv=['mpirun', '-np', f'{self.mpi}', self.abacus])
        calc = Abacus(profile=profile, directory=self.directory,
                      **self.parameters)

        return calc

    def set_neb_chain(self, fmax=0.05, climb=True, interpolate='idpp', traj='neb.traj',):
        """Set neb_chain, namely:
        1. images defining path from initial to final state
        2. neb object including the images and the implemented NEB method
        
        fmax (float): threshold (unit: eV/Angstrom) of the force convergence
        climb (bool): climbing image NEB method
        interpolate (string): interpolate chain path, 'linear' or 'idpp' or None
        """
        # set images
        # continuative running setting
        if Path(traj).exists():
            print(f"----- traj file {traj} exists !!! -----")
            print(f"----- Read Trajectory from {traj}  -----")
            images = read(f"{traj}@-{self.n_max + 2}:")
            interpolate = None
        elif self.guess:
            images = self.guess
            interpolate = None
        else:
            images = [self.initial]
            for i in range(self.n_max):
                image = self.initial.copy()
                images.append(image)
            images.append(self.final)
        # set calculator
        for num, image in enumerate(images):
            # for parallel
            if self.parallel:
                if world.rank == num & world.size:
                    image.calc = self.set_calculator()
            else:
                image.calc = self.set_calculator()

        if self.algorism in ["aseneb", "improvedtangent", "eb", "spline", "string"]:
            if self.dyneb:
                # dynamic neb can only be performed by serial
                print("----- Running Dynamic NEB -----")
                print(f"----- {self.algorism} method is being used -----")
                print("----- Default scale_fmax = 1.0 -----")
                neb = DyNEB(images, climb=climb, fmax=fmax, dynamic_relaxation=True, allow_shared_calculator=True,
                method=self.algorism, parallel=self.parallel, scale_fmax=1.0)
            else:
                print("----- Running ASE-NEB -----")
                print(f"----- {self.algorism} method is being used -----")
                neb = NEB(images, climb=climb, fmax=fmax, method=self.algorism, parallel=self.parallel)
        else:
            print("Error: NEB algorism not supported")
            print("Please choose algorism from 'aseneb', 'improvedtangent', 'eb', 'spline', 'string', ")
            print("AutoNEB method should be used by AbacusAutoNEB")
            exit(1)
        # set interpolate if not set to None
        if interpolate in ['idpp','linear']:
            neb.interpolate(method=interpolate)
            # print-out guess information
            write(f'{interpolate}-guess-{traj}', images, format='traj')
        elif interpolate:
            print("---- Warning: interpolate method not supported, using default linear interpolate ----")
            neb.interpolate(method="linear") # using default
            write(f'linear-guess-{traj}', images, format='traj')
        else:
            print("---- images read from traj file or guess, no interpolate method is used ----")
        return neb

    def run(self, optimizer=FIRE, fmax=0.05, climb=True, interpolate='idpp', traj="neb.traj"):
        """Run Abacus NEB

        optimizer (Optimizer object): defaults to FIRE. BFGS, LBFGS, GPMin, MDMin and QuasiNewton are supported, recommend FIRE method
        fmax (float): threshold (unit: eV/Angstrom) of the force convergence
        climb (bool): climbing image NEB method
        interpolate (string): interpolate chain path, 'linear' or 'idpp' or None, default is 'idpp' for optimal guess chain
        """
        # added continuative running
        neb = self.set_neb_chain(fmax, climb, interpolate, traj)
        if Path(traj).exists():
            traj_restart = f"restart_{traj}"
            print(f"----- OUT Trajectory change to {traj_restart}  -----")
            opt = optimizer(neb, trajectory=traj_restart)
        else:
            opt = optimizer(neb, trajectory=traj)
        opt.run(fmax)

    def _nebtools(self, images):
        return NEBTools(images)

    def get_barriers(self):
        """Returns the barrier estimate from the NEB, along with the Delta E of the elementary reaction"""
        images = read(f'neb.traj@-{self.n_max + 2}:')
        for atoms in images:
            print(atoms.get_potential_energy())

        return self._nebtools(images).get_barrier()

    def plot_bands(self):
        """Given a trajectory containing many steps of a NEB, makes plots of each band in the series in a single PDF"""
        images = read(f'neb.traj@-{self.n_max + 2}:')
        return self._nebtools(images).plot_bands()

# end class AbacusNEB


# start subclass AbacusAutoNEB
# not finished
# autoneb seems to be indenpendent in neb system now
class AbacusAutoNEB(AbacusNEB):
    """Abacus NEB with AutoNEB method"""
    
    def __init__(self, initial, final, parameters, abacus='abacus', algorism="eb", directory='OUT', mpi=1, omp=1, guess=[], n_max=8, parallel=True) -> None:
        super(AbacusAutoNEB, self).__init__(initial, final, parameters, abacus, algorism, directory, mpi, omp, guess, n_max, parallel)
    
    def set_neb_chain(self, fmax=0.05, climb=True, interpolate='idpp'):
        """Set neb_chain, namely:
        1. images defining path from initial to final state
        2. neb object including the images and the implemented NEB method
        
        fmax (float): threshold (unit: eV/Angstrom) of the force convergence
        climb (bool): climbing image NEB method
        interpolate (string): interpolate path linearly from initial to final state, also interpolate can be set to 'idpp'
        """
        if self.guess:
            images = self.guess
            for image in images:
                image.calc = self.set_calculator()
        else:
            images = [self.initial]
            for i in range(self.n_max):
                image = self.initial.copy()
                image.calc = self.set_calculator()
                images.append(image)
            images.append(self.final)
        # need to be specified
        neb = AutoNEB()
        if interpolate:
            neb.interpolate()
        return neb


