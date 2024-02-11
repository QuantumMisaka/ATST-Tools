# using neb chain to get the main vibrate atoms
# JamesMisaka in 2024-02-11

from ase.io import read
from ase.mep.neb import NEBTools
import numpy as np
import sys


def neb2vib(neb_traj):
    """Vibrational analysis from finite displacement by using abacus"""
    neb_traj = read('neb_latest.traj', index=':')
    norm = np.linalg.norm

    # using displacement vector to get the main vibrate atoms
    n_images = n_images = NEBTools(neb_traj)._guess_nimages()
    neb_chain = neb_traj[-n_images:]
    # get TS information from NEB chain
    barrier = NEBTools(neb_chain).get_barrier()[0]
    fmax = NEBTools(neb_chain).get_fmax()
    raw_barrier = max([image.get_potential_energy() for image in neb_chain])
    TS_info = [(ind, image) 
            for ind, image in enumerate(neb_chain) 
            if image.get_potential_energy() == raw_barrier][0]
    step_before_TS = 1
    step_after_TS = 1
    ind_before_TS = TS_info[0] - step_before_TS
    ind_after_TS = TS_info[0] + step_after_TS
    img_TS = TS_info[1]
    img_before = neb_chain[ind_before_TS]
    img_after = neb_chain[ind_after_TS]
    image_vector = (img_before.positions - img_after.positions)
    len_vector = norm(image_vector)
    norm_vector = norm(image_vector / len_vector, axis=1)
    print(norm_vector)
    # those atoms with large contribution to displacement will be picked
    vib_indices = [atom.index for atom in img_TS if norm(norm_vector[atom.index]) > 0.10]
    print(f"=== Locate TS in {TS_info[0]} of 0-{n_images-1} images  ===")
    print(f"=== TS Barrier: {barrier:.4f} (eV) ===")
    print(f"=== TS main moving atoms: {vib_indices} ===")
    return img_TS, vib_indices

if __name__ == "__main__":
    neb_traj = read(sys.argv[1])
    neb2vib(neb_traj)