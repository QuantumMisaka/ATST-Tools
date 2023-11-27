# JamesMisaka in 2023-11-16
# read neb (or dimer) traj file, transform it to 
# 1. extended xyz file, which is widely used in ASE, Ovito, VMD, Nequip and Others, in there will be energy, force and stress information
# 2. ABACUS STRU file, for continuation 
# 3. cif file, for M$ and VESTA to read and analysis
# part of ATST-Tools scripts

import os, sys
from ase.io import read, write, Trajectory
from ase.mep.neb import NEBTools

def write_traj(tag, traj, format="extxyz"):
    """output trajectory file by format"""
    
    if format == "extxyz":
        write(f"{tag}.extxyz", traj, format="extxyz")
    elif format == "traj":
        if os.path.exists(f"{tag}.traj"):
           raise FileExistsError(f"{tag}.traj already exists, You CANNOT transfer a traj file to itself !")
        write(f"{tag}.traj", traj, format="traj")
    elif format ==  "abacus":
        os.mkdir(f"{tag}_STRUs")
        os.chdir(f"{tag}_STRUs")
        for ind, image in enumerate(traj):
            write(f"STRU_{ind:05d}", image, format="abacus")
    elif format == "cif":
        os.mkdir(f"{tag}_CIFs")
        os.chdir(f"{tag}_CIFs")
        for ind, image in enumerate(traj):
            write(f"{tag}-{ind:05d}.cif", image, format="cif")
    else:
        raise ValueError("output format not supported")

if __name__ == "__main__":
    msg = '''
Usage:
    python traj_transform.py [traj_file] [output_format] [optional]
Args:
    traj_file: the trajectory file to be transformed
    output_format: the format of output file, support extxyz, abacus, cif, traj Now. Notice that traj format is only used for --neb option
Optional:
    --neb: to analysis neb traj by neb method, and cut it by band
    '''
    if len(sys.argv) < 3:
        print(msg)
    elif "--neb" in sys.argv[3:]:
        # analysis neb traj by neb method
        print("===> Transform NEB Traj by NEBTools <===")
        traj_file = sys.argv[1]
        tag = traj_file.split(".")[0]
        traj = read(sys.argv[1], ":")
        nebtools = NEBTools(traj)
        image_num = nebtools._guess_nimages()
        band_num = len(traj) // image_num
        output_format = sys.argv[2]
        # cut neb traj by image number and print out
        NEB_dir = "Traj_NEB"
        if not os.path.exists(NEB_dir):
            os.mkdir(NEB_dir)    
        os.chdir(NEB_dir)
        for i in range(band_num):
            neb_band_tag = f"{tag}-{i:04d}"
            band_traj = traj[i*image_num:(i+1)*image_num]
            write_traj(neb_band_tag, band_traj, format=output_format)
        print(f"===> Successfully Transform NEB Traj to {output_format} files by NEB method ! <===")
    else:
        traj_file = sys.argv[1]
        tag = traj_file.split(".")[0]
        traj = read(sys.argv[1], ":")
        output_format = sys.argv[2]
        write_traj(tag, traj, format=output_format)
        print(f"===> Successfully Transform NEB Traj to {output_format} files ! <===")
    
        