#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=32
#SBATCH -J FeC-C-NEB
#SBATCH -o run_neb.out
#SBATCH -e run_neb.err
#SBATCH -p C064M0256G

# JamesMisaka in 2023-11-02
# workflow of ase-abacus-neb method
# each python script do their single job
# ntasks-per-node is for images, cpus-per-task is for calculator
# for one calculator, ntasks-per-node for mpi, cpus-per-task for openmp

# n_image selection should consider the computer resource

# in developer's PKU-WM2 server
source /lustre/home/2201110432/apps/miniconda3/etc/profile.d/conda.sh
conda activate gpaw-intel
module load abacus/3.4.1-icx

# variable
INIT="IS/OUT.ABACUS/running_relax.log"   # running_relax/scf.log
FINAL="FS/OUT.ABACUS/running_relax.log"
#N_IMG=8 # image number in intermediate, each process do one image calculation
# in autoneb, N_IMG is the parallelled-run img number
N_IMG=$SLURM_NTASKS
echo "N_IMG is ${N_IMG}" 

# Job state 
echo $SLURM_JOB_ID > JobRun.state
echo "Start at $(date)" >> JobRun.state

# Job Starting
echo "===== NEB Job Starting at $(date) =====" 

# prepare neb image chain
echo " ===== Make Initial NEB Guess ====="
python neb_make.py $INIT $FINAL $N_IMG --fix 0.14:1 --mag Fe:2.0,C:-0.3 

# neb_make Usage: 
# Default: 
#     python neb_make.py [init_image] [final_image] [n_max] [optional]
# [optional]:
#     --fix [height]:[direction] : fix atom below height (fractional) in direction (0,1,2 for x,y,z)
#     --mag [element1]:[magmom1],[element2]:[magmom2],... : set initial magmom for atoms of element
# Use existing guess: 
#     python neb_make.py -i [init_guess_traj]

# run neb
echo "===== Running NEB ====="

#python neb_run.py
# mpirun -np $N_IMG gpaw python neb_run.py #2>run_neb.err | tee run_neb.out
# if trans-nodes
srun hostname -s | sort -n > slurm.hosts
mpirun -np $N_IMG -machinefile slurm.hosts gpaw python neb_run.py

# post_precessing
echo "===== NEB Process Done ! ====="
echo "===== Running Post-Processing ====="

python neb_post.py neb.traj $N_IMG

echo "===== Done at $(date)! ====="

# Job State
echo "End at $(date)" >> JobRun.state
