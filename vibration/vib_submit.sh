#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=4
#SBATCH -J vib_test
#SBATCH -o running_vib.out
#SBATCH -e running_vib.err
#SBATCH -p C064M0256G
#SBATCH --qos=low

# JamesMisaka in 2023-11-30
# Vibrational analysis by using abacus
# part of ATST-Tools scripts

# in developer's PKU-WM2 server
source /lustre/home/2201110432/apps/miniconda3/etc/profile.d/conda.sh
conda activate ase
module load abacus/3.7.0-icx

echo "Vibrational Calculation Start at $(date) !"
# Job state 
echo $SLURM_JOB_ID > JobRun.state
echo "Start at $(date)" >> JobRun.state

python vib_analysis.py #2>running_vib.err | tee running_vib.out

# post_precessing
echo "===== Done at $(date)! ====="

# Job State
echo "End at $(date)" >> JobRun.state

