#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=4
#SBATCH -J DIMER-ABACUS
#SBATCH -o running_dimer.out
#SBATCH -e running_dimer.err
#SBATCH -p C064M0256G
#SBATCH --qos=normal

# JamesMisaka in 2024-09-16
# workflow of ase-abacus-dimer method

# in developer's PKU-WM2 server
source /lustre/home/2201110432/apps/miniconda3/bin/activate
conda activate ase
module load abacus/3.7.5-icx

# if one just done neb calculation and done neb_post.py
# python neb2dimer.py neb_latest.traj # or neb.traj
# python neb2dimer.py STRU_IS STRU_FS

# Job state 
echo $SLURM_JOB_ID > JobRun.state
echo "Start at $(date)" >> JobRun.state

# just do dimer by run !
# do check your ntasks(mpi) and cpu(omp) setting is right
python neb2dimer_abacus.py

echo "===== NEB2Dimer Calculation Done at $(date)! ====="

# Job State
echo "End at $(date)" >> JobRun.state
