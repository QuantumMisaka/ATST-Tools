# JamesMisaka in 2023-11-02
# Analyze NEB calculation result 
# part of ASE-NEB-ABACUS scripts

from ase.mep.neb import NEB, DyNEB, NEBTools # newest ase
from ase.io import read, write
import sys, os

# next dev: point-out transition state
class NEBPost():
    """Post-Processing of NEB result from traj file"""
    def __init__(self, traj_file, n_max):
        self.all_image = read(traj_file, index=":", format="traj")
        self.n_images = n_max + 2
        self.neb_chain = read(traj_file, index=f"{- self.n_images}:", format="traj")

    def get_barrier(self):
        """Returns energy of all image and the barrier estimate from the NEB, along with the Delta E of the elementary reaction"""
        for i, atoms in enumerate(self.neb_chain):
            ene = atoms.get_potential_energy()
            print(f"num: {i}; Energy: {ene} (eV)")
        barrier = NEBTools(self.neb_chain).get_barrier()
        print(f"Reaction Barrier and Energy Difference: {barrier} (eV)")
        return barrier

    def plot_neb_bands(self, label='nebplots_chain'):
        """makes plots of final neb band in the series in a single PDF"""
        return NEBTools(self.neb_chain).plot_bands(nimages=self.n_images, label=label)
    
    def plot_all_bands(self, label='nebplots_all'):
        """Gmakes plots of all band during neb in the series in a single PDF"""
        return NEBTools(self.all_image).plot_bands(label=label)
    
    def view_neb_bands(self):
        """view neb chain"""
        os.system(f'ase gui neb.traj@-{self.n_images}:')

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python neb_postprocess.py [traj_file] [n_images]")
    else:
        traj_file = sys.argv[1]
        n_max = int(sys.argv[2])
        result = NEBPost(traj_file, n_max)
        result.get_barrier()
        result.plot_all_bands()
        result.plot_neb_bands()
    