# JamesMisaka in 2023-12-04
# read STRU or other atoms files, transfer them to a traj file 
# for doing continuation calculation, espeacially for AutoNEB
# part of ATST-Tools scripts

import sys
from ase.io import read, write
from ase.atoms import Atoms

def traj_collect(files: list[Atoms], out_file="collection.traj", no_calc=False):
    """Collect Atoms files to a traj file contain list[Atoms]
    
    Args:
        files (list): list of Atoms contrain all atom files
        out_file (string): output file name
        no_calc (bool): if output trajectory have property and calculator, default no
    """
    traj = []
    if no_calc:
        traj = [atoms.copy() for atoms in files]
    else:
        traj = files
    write(out_file, traj, format="traj")
    print("==> Successfully Collect All STRU files to Trajectory file <==")
    
if __name__ == "__main__":
    msg = '''
Usage: 
    python traj_collect.py  [(--no-calc)]  [structures]
structures list should contrain 2 or more structures
    --no-calc is optional for not keeping properties
'''
    if len(sys.argv) <= 2:
        print(msg)
    else:
        no_calc = False
        if "--no-calc" in sys.argv:
            if len(sys.argv) <= 3:
                print(msg)
                exit()
            index_no_calc = sys.argv.index("--no-calc")
            no_calc = True
            stru_files = sys.argv[1:index_no_calc]
            stru_files.extend(sys.argv[index_no_calc+1:])
        else:
            stru_files = sys.argv[1:]
        try:
            atoms_list = [read(stru) for stru in stru_files]
        except:
            raise ValueError("Can't read structures files")
        traj_collect(atoms_list, no_calc=no_calc)
        print("==> Done ! <==")
