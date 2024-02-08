import numpy as np
import os

from ase.io import read, write, Trajectory
from ase import Atoms
from ase.optimize import BFGS, FIRE, QuasiNewton
from ase.constraints import FixAtoms
from ase.visualize import view
from ase.mep.neb import NEBTools, NEB, DyNEB
from ase.mep.autoneb import AutoNEB
from ase.mep.dimer import DimerControl, MinModeAtoms, MinModeTranslate
from deepmd_pt.utils.ase_calc import DPCalculator as DP


from ase.mep.autoneb import AutoNEB

model = "FeCHO-dpa2-full.pt"
n_max = 8
neb_fmax = 1.00  # neb should be rough
dimer_fmax = 0.05 # dimer use neb guess
climb = True
scale_fmax = 1.0 # use dyneb to reduce message far from TS
omp = 16
neb_algorism = "improvedtangent"
neb_traj = "neb_images.traj"

os.environ['OMP_NUM_THREADS'] = "omp"

# init and final stru
atom_init = read("./C-C-Fe5C2_510/data/IS/STRU")
atom_final = read("./C-C-Fe5C2_510/data/FS/STRU")
atom_init.calc = DP(model=model)
atom_final.calc = DP(model=model)
init_relax = BFGS(atom_init)
final_relax = BFGS(atom_final)
init_relax.run(fmax=0.05)
final_relax.run(fmax=0.05)

write("init_opted.traj", atom_init, format="traj")
write("final_opted.traj", atom_final, format="traj")

# run neb and dimer 
# function setting

class DPDimer:
    """Customize Dimer calculation workflow by using Deep Potential"""
    
    def __init__(self, init_Atoms, model,
                 omp=1, directory='DIMER', 
                 traj_file='dimer.traj',
                 init_eigenmode_method='displacement',
                 displacement_vector: np.ndarray = None,):
        """Initialize Dimer method by using ASE-ABACUS

        init_Atoms (Atoms object): starting image, can be from every way including NEB result
        parameters (dict): settings of abacus input parameters
        model (str): DeepPotential model
        directory (str): calculator directory name, for parallel calculation {directory}-rank{i} will be the directory name
        mpi (int): number of MPI for abacus calculator
        omp (int): number of OpenMP for abacus calculator
        traj_file (str): trajectory file name for dimer calculation, when running dimer calculation, trajetory will be written to this file, default is 'dimer.traj'
        init_eigenmode_method (str): dimer initial eigenmode method. Choose from 'displacement' and 'gauss'.
        displacement_vector (np.ndarray): displacement vector for dimer initial eigenmode. Only used when init_eigenmode_method is 'displacement'
        """
        self.init_Atoms = init_Atoms
        self.model = model
        self.omp = omp
        self.directory = directory
        self.traj_file = traj_file
        self.init_eigenmode_method = init_eigenmode_method
        self.displacement_vector = displacement_vector
        
    def set_calculator(self):
        """Set Abacus calculators"""
        os.environ['OMP_NUM_THREADS'] = f'{self.omp}'
        calc = DP(model=self.model)
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
        d_mask = [True] * len(dimer_init)
        const = dimer_init._get_constraints()
        # const will be empty list if no constraint
        if const:
            const_object = dimer_init._get_constraints()[0].get_indices()
            for ind in const_object:
                d_mask[ind] = False
            return d_mask
        else:
            print("--- Notice: No constraint found in init Atoms, there will be no mask in dimer calculation ---")
            return d_mask
    
    def set_d_mask_by_specified(self, moving_atoms_ind: list):
        """set mask be choosing moving atoms, the others are masked"""
        print(f"=== Set mask by specifing moving atoms {moving_atoms_ind} ===")
        dimer_init = self.init_Atoms
        d_mask = [False] * len(dimer_init)
        for ind in moving_atoms_ind:
            d_mask[ind] = True
        return d_mask
        
    def run(self, fmax=0.05, properties=["energy", "forces", "stress"], moving_atoms_ind: list = None):
        """run dimer calculation workflow
        
        Args:
            fmax (float): threshold (unit: eV/Angstrom) of the force convergence
            properties (list): properties dumped in trajectory files, default ['energy', 'forces', 'stress']
        """
        dimer_init = self.init_Atoms
        dimer_init.calc = self.set_calculator()
        dimer_traj = Trajectory(self.traj_file, 'w', dimer_init, properties=properties)
        if self.init_eigenmode_method == "displacement":
            # d_mask = self.set_d_mask_by_displacement() # not the best
            if moving_atoms_ind:
                d_mask = self.set_d_mask_by_specified(moving_atoms_ind)
            else:
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


# run neb
images = [atom_init]
for i in range(n_max):
    image = atom_init.copy()
    image.set_calculator(DP(model=model))
    images.append(image)
images.append(atom_final)
neb = DyNEB(images, 
            climb=climb, dynamic_relaxation=True, fmax=neb_fmax,
            method=neb_algorism, parallel=False, scale_fmax=scale_fmax)
neb.interpolate(method="idpp")

traj = Trajectory(neb_traj, 'w', neb)
opt = FIRE(neb, trajectory=traj)
opt.run(neb_fmax)

# neb displacement to dimer
n_images = NEBTools(images)._guess_nimages()
neb_raw_barrier = max([image.get_potential_energy() for image in images])
fmax = NEBTools(images).get_fmax()
barrier = NEBTools(images).get_barrier()[0]
TS_info = [(ind, image) 
            for ind, image in enumerate(images) 
            if image.get_potential_energy() == neb_raw_barrier][0]
print(f"=== Locate TS in {TS_info[0]} of 0-{n_images-1} images  ===")
print(f"=== NEB Raw Barrier: {neb_raw_barrier:.4f} (eV) ===")
print(f"=== NEB Fmax: {fmax:.4f} (eV/A) ===")
print(f"=== Now Turn to Dimer with NEB Information ===")

# para for neb2dimer
step_before_TS = 1
step_after_TS = 1
norm_vector = 0.01
#out_vec = 'displacement_vector.npy',

ind_before_TS = TS_info[0] - step_before_TS
ind_after_TS = TS_info[0] + step_after_TS
img_before = images[ind_before_TS]
img_after = images[ind_after_TS]
image_vector = (img_before.positions - img_after.positions)
modulo_norm = np.linalg.norm(image_vector) / norm_vector
displacement_vector = image_vector / modulo_norm
print(f"=== Displacement vector generated by {ind_before_TS} and {ind_after_TS} images of NEB chain ===")
print(f"=== Which is normalized to {norm_vector} length ! ===")
#np.save(out_vec,displacement_vector)

# dimer part
dimer_init = TS_info[1].copy()
init_eigenmode_method = "displacement"
dimer = DPDimer(dimer_init, model=model,
                        omp=omp, 
                        init_eigenmode_method=init_eigenmode_method,
                        displacement_vector=displacement_vector)
dimer.run(fmax=dimer_fmax)