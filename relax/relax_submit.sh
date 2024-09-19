#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=4
#SBATCH -J RELAX-ASE-ABACUS
#SBATCH -o running_relax.out
#SBATCH -e running_relax.err
#SBATCH -p C064M0256G
#SBATCH --qos=normal

# JamesMisaka in 2023-11-14
# workflow of ase-abacus-dimer method
# for one calculator, ntasks-per-node for mpi, cpus-per-task for openmp

# in developer's PKU-WM2 server
source /data/softwares/miniconda3/bin/activate
conda activate gpaw-intel
module load abacus/3.7.0-icx

# Job state 
echo $SLURM_JOB_ID > JobRun.state
echo "Start at $(date)" >> JobRun.state

# just do dimer by run !
# do check your ntasks(mpi) and cpu(omp) setting is right
python relax_run.py

echo "===== Relax Calculation Done at $(date)! ====="

# Job State
echo "End at $(date)" >> JobRun.state
