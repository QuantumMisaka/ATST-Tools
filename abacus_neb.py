import os 
from ase.neb import NEB, DyNEB
from ase.autoneb import AutoNEB
from ase.calculators.abacus import Abacus, AbacusProfile
from ase.optimize import FIRE, BFGS
from ase.neb import NEBTools
from ase.io import read, write
from pathlib import Path

# os.environ['ABACUS_PP_PATH'] = '/home/james/example/PP'
# os.environ['ABACUS_ORBITAL_PATH'] = '/home/james/example/ORB'
# print(os.environ.get('ABACUS_PP_PATH'))
# print(os.environ.get('ABACUS_ORBITAL_PATH'))


class AbacusNEB:
    """Customize Abacus workflow for Nudged Elastic Band calculation"""

    def __init__(self, initial, final, parameters, abacus='abacus', algorism="dyneb", directory='OUT', mpi=1, omp=1, guess=[], n_max=8) -> None:
        """Initialize initial and final states

        initial (Atoms object): relaxed starting image
        final (Atoms object): relaxed ending image
        n_max (int): the number of images to interpolate between the initial and final images.
        parameters (dict): settings of input parameters
        abacus (str): Abacus executable file. Default: 'abacus'
        algorism (str): NEB algorism. which can be
        - 'aseneb'
        - 'improvedtangent'
        - 'eb'
        - 'spline'
        - 'string'
        - 'autoneb'
        - 'dyneb'
        Default: 'dyneb'
        
        directory (str): work directory
        mpi (int): number of MPI
        omp (int): number of OpenMP
        guess (List(Atoms)): intermediate guesses. If set, not need to set n_max
        """

        self.initial = initial
        self.final = final
        self.n_max = n_max
        self.parameters = parameters
        self.algorism = algorism
        self.abacus = abacus
        self.directory = directory
        self.mpi = mpi
        self.omp = omp
        self.guess = guess

    def set_calculator(self):
        """Set Abacus calculators"""
        os.environ['OMP_NUM_THREADS'] = f'{self.omp}'
        profile = AbacusProfile(
            argv=['mpirun', '-np', f'{self.mpi}', self.abacus])
        calc = Abacus(profile=profile, directory=self.directory,
                      **self.parameters)

        return calc

    def set_neb_chain(self, fmax=0.05, climb=True, interpolate=True):
        """Set neb_chain, namely:
        1. images defining path from initial to final state
        2. neb object including the images and the implemented NEB method
        
        fmax (float): threshold (unit: eV/Angstrom) of the force convergence
        climb (bool): climbing image NEB method
        interpolate (bool): interpolate path linearly from initial to final state
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
        if self.algorism == "dyneb":
            neb = DyNEB(images, climb=climb, fmax=fmax, dynamic_relaxation=True, allow_shared_calculator=True)
        elif self.algorism in ["aseneb", "improvedtangent", "eb", "spline", "string"]:
            neb = NEB(images, climb=climb, fmax=fmax, method=self.algorism, allow_shared_calculator=True)
        else:
            print("Error: NEB algorism not supported")
            print("Please choose algorism from 'dyneb', 'aseneb', 'improvedtangent', 'eb', 'spline', 'string', 'autoneb', ")
            print("AutoNEB method should be used by AbacusAutoNEB")
            exit(1)
        if interpolate:
            neb.interpolate()
        
        return neb

    def run(self, optimizer=FIRE, fmax=0.05, climb=True, interpolate=True):
        """Run Abacus NEB

        optimizer (Optimizer object): defaults to FIRE. BFGS, LBFGS, GPMin, MDMin and QuasiNewton are supported, recommend FIRE method
        fmax (float): threshold (unit: eV/Angstrom) of the force convergence
        climb (bool): climbing image NEB method
        interpolate (bool): interpolate path linearly from initial to final state
        """
        neb = self.set_neb_chain(fmax, climb, interpolate)
        opt = optimizer(neb, trajectory='neb.traj', logfile='opt.log')
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


if __name__ == '__main__':
    # Run a example: Au diffusion on Al(100) surface
    n_max = 5
    mpi = 64
    omp = 1
    abacus = 'abacus'
    pseudo_dir = "/data/home/liuzq/example/PP"
    basis_dir = "/data/home/liuzq/example/ORB"
    pp = {"Al": "Al_ONCV_PBE-1.0.upf",
          "Au": "Au_ONCV_PBE-1.0.upf", }
    basis = {"Al": "Al_gga_7au_100Ry_4s4p1d.orb",
             "Au": "Au_gga_7au_100Ry_4s2p2d1f.orb"}
    kpts = [2, 2, 1]
    parameters = {
        'calculation': 'scf',
        'xc': 'pbe',
        'ecutwfc': 100,
        'smearing_method': 'gaussian',
        'smearing_sigma': 0.01,
        'basis_type': 'lcao',
        'ks_solver': 'genelpa',
        'mixing_type': 'pulay',
        'scf_thr': 1e-8,
        'out_chg': 1,
        'out_bandgap': 1,
        'kpts': kpts,
        'pp': pp,
        'basis': basis,
        'pseudo_dir': pseudo_dir,
        'basis_dir': basis_dir,
        'vdw_method': 'd3_bj',
        'cal_force': 1,
    }

    from ase.build import fcc100, add_adsorbate
    from ase.constraints import FixAtoms
    from ase.optimize import QuasiNewton
    slab = fcc100('Al', size=(2, 2, 3))
    add_adsorbate(slab, 'Au', 1.7, 'hollow')
    slab.center(axis=2, vacuum=4.0)
    mask = [atom.tag > 1 for atom in slab]
    slab.set_constraint(FixAtoms(mask=mask))
    os.environ['OMP_NUM_THREADS'] = f'{omp}'
    profile = AbacusProfile(
        argv=['mpirun', '-np', f'{mpi}', abacus])
    slab.calc = Abacus(profile=profile, directory="INIT",
                      **parameters)
    
    qn = QuasiNewton(slab, trajectory='initial.traj')
    qn.run(fmax=0.05)
    slab[-1].x += slab.get_cell()[0, 0] / 2
    qn = QuasiNewton(slab, trajectory='final.traj')
    qn.run(fmax=0.05)
    # Initial state:
    initial = read('initial.traj')

    # Final state:
    final = read('final.traj')
    
    neb = AbacusNEB(initial=initial, final=final, parameters=parameters,
                   mpi=mpi, omp=omp, abacus=abacus, n_max=n_max)
    neb.run(optimizer=BFGS, climb=False, interpolate=True, fmax=0.05)

    # Get barrier
    barrier = neb.get_barriers()
    print(barrier)
    neb.plot_bands()

    # Visualize the results
    # os.system(f'ase gui neb.traj@-{n_max}:')

