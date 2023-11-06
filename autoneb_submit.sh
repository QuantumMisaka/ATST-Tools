#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH -J AutoNEB-ABACUS
#SBATCH -o run_autoneb.out
#SBATCH -e run_autoneb.err
#SBATCH -p C064M0256G

# JamesMisaka in 2023-11-02
# workflow of abacus-autoneb method
# each python script do their single job
# ntasks-per-node is for n_simul, cpus-per-task is for calculator
# for one calculator, ntasks-per-node for mpi, cpus-per-task for openmp

# in developer's PKU-WM2 server
source /lustre/home/2201110432/apps/miniconda3/etc/profile.d/conda.sh
conda activate gpaw-intel
module load abacus/3.4.2-icx

# variable
INIT="INIT/OUT.ABACUS/running*.log"   # running_relax/scf.log
FINAL="FINAL/OUT.ABACUS/running*.log"
#NSIMUL=4 # the parallelled-run img number
NSIMUL=$SLURM_NTASKS
echo "NSIMUL is ${NSIMUL}" 

# Job state 
echo $SLURM_JOB_ID > JobRun.state
echo "Start at $(date)" >> JobRun.state

# Job Starting
echo "===== AutoNEB Job Starting =====" 

# prepare neb image chain
echo " ===== Make Initial NEB Guess ====="
python neb_make.py $INIT $FINAL 0 # --fix 0.2:1 --mag Fe:1.7,C:-0.3 

# neb_make Usage: 
# Default: 
#     python neb_make.py [init_image] [final_image] [n_max] [optional]
# [optional]:
#     --fix [height]:[direction] : fix atom below height (fractional) in direction (0,1,2 for x,y,z)
#     --mag [element1]:[magmom1],[element2]:[magmom2],... : set initial magmom for atoms of element
# Use existing guess: 
#     python neb_make.py -i [init_guess_traj]

# run neb
echo "===== Running AutoNEB ====="

#python neb_run.py
mpirun -np $NSIMUL gpaw python autoneb_run.py #2>run_neb.err | tee run_neb.out
# if trans-nodes
# srun hostname -s | sort -n > slurm.hosts
# mpirun -np $NSIMUL -machinefile slurm.hosts gpaw python autoneb_run.py

# post_precessing
echo "===== AutoNEB Process Done ! ====="
echo "===== Running Post-Processing ====="

python neb_post.py --autoneb run_autoneb???.traj

echo "===== Done ! ====="

# Job State
echo "End at $(date)" >> JobRun.state
