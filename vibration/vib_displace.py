from ase.vibrations import Vibrations
from ase.io import read, write
import os

stru = "STRU_ts"
atoms = read(stru)
# indices setting
vib_indices = list(range(48, 54))
vib_indices += [atom.index for atom in atoms if atom.symbol == 'H']

parent_dir = "disp_test"


pp = {
      'C':'C_ONCV_PBE-1.0.upf',
      'H':'H_ONCV_PBE-1.0.upf',
      'Pt':'Pt_ONCV_PBE-1.0.upf',
      }
basis = {
         'C': 'C_gga_7au_100Ry_2s2p1d.orb',
         'H': 'H_gga_6au_100Ry_2s1p.orb',
         'Pt': 'Pt_gga_7au_100Ry_4s2p2d1f.orb'
         ,}

vib = Vibrations(atoms, indices=vib_indices)

if not os.path.isdir(parent_dir):
    os.mkdir(parent_dir)
os.chdir(parent_dir)

for ind, info in enumerate(vib.iterdisplace()):
    name, image = info
    directory = f"{name}-{ind:03d}"
    if not os.path.isdir(directory):
        os.mkdir(directory)
    write(f"{directory}/STRU", image, format="abacus", pp=pp, basis=basis)