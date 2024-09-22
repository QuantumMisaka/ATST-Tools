# read STRU file, run vibrational analysis jobs
import sys
from ase.vibrations import Vibrations
from ase.thermochemistry import HarmonicThermo
from ase.io import read, write
from ase.io import Trajectory
# from deepmd_pt.utils.ase_calc import DPCalculator as DP
from deepmd.calculator import DP

# setting
model = "FeCHO-dpa2-full.pt"
vib_name = "vib"
vib_indices = None
T = 523.15 # K

# developer only
delta = 0.01
nfree = 2

if len(sys.argv) > 1:
    stru = read(sys.argv[1], format='abacus')
else:
    stru = read("STRU", format='abacus')
if "--ind" in sys.argv[2:]:
    fix_ind = sys.argv.index("--ind")
    vib_indices = [int(i) for i in sys.argv[fix_ind+1:]]


if __name__ == "__main__":
    print("==> Starting Vibrational Analysis <==")

    stru.calc = DP(model=model)

    vib = Vibrations(stru, indices=vib_indices, 
                    name=vib_name, delta=delta, nfree=nfree)

    print("==> Running Vibrational Analysis <==")
    vib.run()
    # post-processing
    print("==> Done !!! <==")
    print(f"==> All force cache will be in {vib_name} directory <==")
    print("==> Vibrational Analysis Summary <==")
    vib.summary()
    print("==> Writing All Mode Trajectory <==")
    vib.write_mode()
    # thermochemistry
    print("==> Doing Harmonmic Thermodynamic Analysis <==")
    vib_energies = vib.get_energies()
    #real_vib_energies = np.array([energy for energy in vib.get_energies() if energy.imag == 0 and energy.real > 0], dtype=float)
    thermo = HarmonicThermo(vib_energies, ignore_imag_modes=True,)
    entropy = thermo.get_entropy(T)
    free_energy = thermo.get_helmholtz_energy(T)
    print(f"==> Entropy: {entropy:.6e} eV/K <==")
    print(f"==> Free Energy: {free_energy:.6f} eV <==")
    print()