#!/bin/bash
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH -J NEB-PtCO-par
#SBATCH -o run_neb.out
#SBATCH -e run_neb.err
#SBATCH -p C064M0256G

# JamesMisaka in 2023-11-02
# workflow of ase-abacus-neb method
# each python script do their single job
# ntasks-per-node is for images, cpus-per-task is for calculator

source /lustre/home/2201110432/apps/miniconda3/etc/profile.d/conda.sh
conda activate gpaw-intel
module load abacus/3.4.1-dev-icx

# variable
INIT="INIT/OUT.ABACUS/running*.log"   # running_scf/relax.log
FINAL="FINAL/OUT.ABACUS/running*.log"
NMAX=8 # image number in intermediate, each process do one image calculation
#NMAX= $SLURM_NTASKS
echo "NMAX is ${NMAX}" 

# Job state 
echo $SLURM_JOB_ID > JobRun.state
echo "Start at $(date)" >> JobRun.state

# prepare neb image chain
echo " ===== Make Initial NEB Guess ====="
python neb_make.py $INIT $FINAL $NMAX --fix 0.10:1 

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
mpirun -np $NMAX gpaw python neb_run.py #2>run_neb.err | tee run_neb.out
# if trans-nodes
# srun hostname -s | sort -n > slurm.hosts
# mpirun -np $NMAX -machinefile slurm.hosts gpaw python neb_run.py

# post_precessing
echo "===== NEB Process Done ! ====="
echo "===== Running Post-Processing ====="

python neb_post.py neb.traj $NMAX

echo "===== Done ! ====="

# Job State
echo "End at $(date)" >> JobRun.state
