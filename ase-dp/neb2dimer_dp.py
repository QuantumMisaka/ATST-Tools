# Author: JamesMisaka 
# Using deepmd to search TS via NEB-DIMER
# Last Update: 2024-09-17

import numpy as np
import os, sys, shutil
from copy import deepcopy

from ase.io import read, write, Trajectory
from ase import Atoms
from ase.optimize import BFGS, FIRE, QuasiNewton
from ase.constraints import FixAtoms
from ase.visualize import view
from ase.mep.neb import NEBTools, NEB, DyNEB
from ase.mep.dimer import DimerControl, MinModeAtoms, MinModeTranslate
from ase.vibrations import Vibrations
from ase.thermochemistry import HarmonicThermo

from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.diffusion.neb.pathfinder import IDPPSolver

from deepmd.calculator import DP
#from deepmd_pt.utils.ase_calc import DPCalculator as DP


# parameter setting
model = "FeCHO-dpa220-100-30656.pt"
n_neb_images = 12
neb_fmax = 1.00  # neb should be rough
dimer_fmax = 0.05 # dimer use neb guess
climb = True
scale_fmax = 0.0 # 0 or 1, 1 for using dyneb
omp = 16
neb_algorism = "improvedtangent"
neb_log = "neb_images.traj"
dimer_log = "dimer_images.traj"
OPTSolver = QuasiNewton
NEBSolver = FIRE



# developer only
os.environ['OMP_NUM_THREADS'] = "omp"
neb_sort_tol = 0.5 # set 0 to shut down
step_before_TS = 1
step_after_TS = 1
norm_vector = 0.01

# reading part
msg = '''
Usage: 
- For using IS and FS: 
    python neb2dimer_dpa2.py [init_stru] [final_stru] ([format])
- For using existing NEB: 
    python neb2dimer_dpa2.py [neb_latest.traj]
'''
if len(sys.argv) < 2:
    print(msg)
    sys.exit(1)
elif len(sys.argv) == 2:
    if sys.argv[1] == "-h" or sys.argv[1] == "--help":
        print(msg)
        sys.exit(0)
    else:
        neb_traj_read = sys.argv[1]
        neb_abacus = read(neb_traj_read, ":", format="traj")
        atom_init = neb_abacus[0]
        atom_final = neb_abacus[-1]
        assert type(atom_init) == Atoms and type(atom_final) == Atoms, \
        "The input file is not a trajectory file contained Atoms object"
else:
    init_stru = sys.argv[1]
    final_stru = sys.argv[2]
    if len(sys.argv) == 4:
        format = sys.argv[3]
    else:
        format = None # auto detect
    atom_init = read(init_stru, format=format)
    atom_final = read(final_stru, format=format)

atom_init.calc = DP(model=model)
atom_final.calc = DP(model=model)
init_relax = OPTSolver(atom_init)
final_relax = OPTSolver(atom_final)
init_relax.run(fmax=0.05)
final_relax.run(fmax=0.05)

init_opted = "init_opted.traj"
final_opted = "final_opted.traj"
write(init_opted, atom_init, format="traj")
write(final_opted, atom_final, format="traj")

# run neb and dimer 
# function setting

class DPDimer:
    """Customize Dimer calculation workflow by using Deep Potential"""
    
    def __init__(self, init_Atoms, model,
                 omp=1, directory='DIMER', 
                 traj_file='dimer.traj',
                 init_eigenmode_method='displacement',
                 displacement_vector: np.ndarray = None,):
        """Initialize Dimer method by using ASE-DP

        init_Atoms (Atoms object): starting image, can be from every way including NEB result
        parameters (dict): settings of abacus input parameters
        model (str): DeepPotential model
        directory (str): calculator directory name, for parallel calculation {directory}-rank{i} will be the directory name
        omp (int): number of OpenMP for DP calculator
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
        """set mask by constraint of Atoms
        
        Notice: This function have some problem in dealing with abacus STRU, which FixCatesian object will be independent
        """
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
        
    def run(self, fmax=0.05, properties=["energy", "forces", "stress"]):
        """run dimer calculation workflow
        
        Args:
            fmax (float): threshold (unit: eV/Angstrom) of the force convergence
            properties (list): properties dumped in trajectory files, default ['energy', 'forces', 'stress']
        """
        dimer_init = self.init_Atoms
        dimer_init.calc = self.set_calculator()
        dimer_traj = Trajectory(self.traj_file, 'w', dimer_init, properties=properties)
        if self.init_eigenmode_method == "displacement":
            d_mask = self.set_d_mask_by_constraint()
            # d_mask = self.set_d_mask_by_displacement()
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
        
def main4dis(displacement_vector, thr=0.10):
    """Get Main Parts of Displacement Vector by using threshold"""
    len_vector = np.linalg.norm(displacement_vector)
    norm_vector = np.linalg.norm(displacement_vector / len_vector, axis=1)
    main_indices = [ind for ind,vec in enumerate(norm_vector) if vec > thr]
    return main_indices, norm_vector

def thermo_analysis(atoms, T, name="vib", indices=None, delta=0.01, nfree=2):
    """Do Thermo Analysis by using ASE"""
    vib_mode_dir = f"{name}_mode"
    if os.path.exists(name):
        shutil.move(name, f"{name}_bak")
    if not os.path.exists(vib_mode_dir):
        os.mkdir(vib_mode_dir)   
    vib = Vibrations(atoms, indices=indices, name=name, delta=delta, nfree=nfree)
    vib.run()
    vib.summary()
    ROOT_DIR = os.getcwd()
    os.chdir(vib_mode_dir)
    vib.write_mode()
    os.chdir(ROOT_DIR)
    vib_energies = vib.get_energies()
    thermo = HarmonicThermo(vib_energies, ignore_imag_modes=True,)
    entropy = thermo.get_entropy(T)
    free_energy = thermo.get_helmholtz_energy(T)
    print(f"==> Entropy: {entropy:.6e} eV/K <==")
    print(f"==> Free Energy: {free_energy:.6f} eV <==")
    print()

# run neb
# use pymatgen interpolate, which is better and avoid trans-cell problem
print(f'Generating NEB path, number of ase_path: {n_neb_images}, sort_tol: {neb_sort_tol}')
is_pmg = AseAtomsAdaptor.get_structure(atom_init)
fs_pmg = AseAtomsAdaptor.get_structure(atom_final)
path = IDPPSolver.from_endpoints([is_pmg, fs_pmg], n_neb_images, sort_tol=neb_sort_tol)
new_path = path.run(maxiter=5000, tol=1e-5, gtol=1e-3)
# conver path to ase format, and add SinglePointCalculator
ase_path = [i.to_ase_atoms() for i in new_path]
ase_path[0] = atom_init
ase_path[-1] = atom_final
for img in ase_path[1:-1]:
    img.calc = DP(model=model)

neb = DyNEB(ase_path, 
            climb=climb, dynamic_relaxation=True, fmax=neb_fmax,
            method=neb_algorism, parallel=False, scale_fmax=scale_fmax,
            allow_shared_calculator=True)

traj = Trajectory(neb_log, 'w', neb)
opt = NEBSolver(neb, trajectory=traj)
opt.run(neb_fmax)

# neb displacement to dimer
n_images = NEBTools(ase_path)._guess_nimages()
neb_raw_barrier = max([image.get_potential_energy() for image in ase_path])
fmax = NEBTools(ase_path).get_fmax()
barrier = NEBTools(ase_path).get_barrier()[0]
TS_info = [(ind, image) 
            for ind, image in enumerate(ase_path) 
            if image.get_potential_energy() == neb_raw_barrier][0]
print(f"=== Locate TS in {TS_info[0]} of 0-{n_images-1} ase_path  ===")
print(f"=== NEB Raw Barrier: {neb_raw_barrier:.4f} (eV) ===")
print(f"=== NEB Fmax: {fmax:.4f} (eV/A) ===")
print(f"=== Now Turn to Dimer with NEB Information ===")

# displacement vector for vib analysis
#out_vec = 'displacement_vector.npy',

ind_before_TS = TS_info[0] - step_before_TS
ind_after_TS = TS_info[0] + step_after_TS
img_before = ase_path[ind_before_TS]
img_after = ase_path[ind_after_TS]
image_vector = (img_before.positions - img_after.positions)
modulo_norm = np.linalg.norm(image_vector) / norm_vector
displacement_vector = image_vector / modulo_norm
print(f"=== Displacement vector generated by {ind_before_TS} and {ind_after_TS} ase_path of NEB chain ===")
print(f"=== Which is normalized to {norm_vector} length ! ===")
#np.save(out_vec,displacement_vector)

# dimer part
dimer_init = TS_info[1].copy()
init_eigenmode_method = "displacement"
dimer = DPDimer(dimer_init, model=model,
                        omp=omp, 
                        init_eigenmode_method=init_eigenmode_method,
                        traj_file=dimer_log,
                        displacement_vector=displacement_vector)
dimer.run(fmax=dimer_fmax)

# get struc of IS,FS,TS
write("IS_get.cif", atom_init, format="cif")
write("FS_get.cif", atom_final, format="cif")
write("TS_get.cif", dimer_init, format="cif")
write("IS_get.stru", atom_init, format="abacus")
write("FS_get.stru", atom_final, format="abacus")
write("TS_get.stru", dimer_init, format="abacus")

# get energy informations
ene_init = atom_init.get_potential_energy()
ene_final = atom_final.get_potential_energy()
ene_ts = dimer_init.get_potential_energy()
ene_delta = ene_final - ene_init
ene_activa = ene_ts - ene_init
ene_act_rev = ene_ts - ene_final
msg = f'''
==> TS-Search Results <==
- Items      Energy
- IS         {ene_init:.6f}
- FS         {ene_final:.6f}
- TS         {ene_ts:.6f}
- dE         {ene_delta:.6f}
- Ea_f       {ene_activa:.6f}
- Ea_r       {ene_act_rev:.6f}
'''
print(msg)

# use neb2dimer information to do vibration analysis
print("==> Do Vibrational Analysis by DP Potential <==")
vib_indices, norm_vector = main4dis(image_vector, thr=0.10)
print(f"=== TS main moving atoms: {vib_indices} ===")
T = 523.15 # K
delta = 0.01
nfree = 2

vib_is_name = 'vib_is'
vib_fs_name = 'vib_fs'
vib_ts_name = 'vib_ts'

print("==> For TS Structure <==")
thermo_analysis(dimer_init, T, name=vib_ts_name, indices=vib_indices, delta=delta, nfree=nfree)
print("==> For Initial Structure <==")
thermo_analysis(atom_init, T, name=vib_is_name, indices=vib_indices, delta=delta, nfree=nfree)
print("==> For Final Structure <==")
thermo_analysis(atom_final, T, name=vib_fs_name, indices=vib_indices, delta=delta, nfree=nfree)



