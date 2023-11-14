import numpy as np
import sys
from ase.io import Trajectory, read, write


def neb2dimer(neb_chain, out_traj='dimer_init.traj', out_vec='displacement_vector.npy'):
    # use neb traj for dimer
    # you should use the latest neb chain or other single neb chain
    neb_chain = Trajectory()
    raw_barrier = max([image.get_potential_energy() for image in neb_chain])
    # get TS information from NEB chain
    TS_info = [(ind, image) 
            for ind, image in enumerate(neb_chain) 
            if image.get_potential_energy() == raw_barrier][0]
    # output TS of neb for dimer init
    write(out_traj, TS_info[1], format='traj')
    # output displancement vector
    img_before_TS = neb_chain[TS_info[0] - 1]
    img_after_TS = neb_chain[TS_info[0] + 1]
    displacement_vector = (img_after_TS.positions - img_before_TS.positions)
    np.save(out_vec,displacement_vector)
    return
    
if __name__ == "__main__":
    msg = '''
    
    '''