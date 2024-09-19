#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=32
#SBATCH -J NEB-ABACUS
#SBATCH -o running_neb.out
#SBATCH -e running_neb.err
#SBATCH -p C064M0256G
#SBATCH --qos=normal

# JamesMisaka in 2023-11-02
# workflow of ase-abacus-neb method
# each python script do their single job
# ntasks-per-node is for images, cpus-per-task is for calculator
# for one calculator, ntasks-per-node for mpi, cpus-per-task for openmp

# n_image selection should consider the computer resource

# in developer's PKU-WM2 server
source /lustre/home/2201110432/apps/miniconda3/bin/activate
conda activate gpaw-intel
module load abacus/3.7.0-icx

# variable
INIT="INIT/OUT.ABACUS/running*.log"   # running_relax/scf.log
FINAL="FINAL/OUT.ABACUS/running*.log"
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
# one may do this in login node better
echo " ===== Make Initial NEB Guess ====="
python neb_make.py -n $NSIMUL -i $INIT $FINAL # --fix 0.2:1 --mag Fe:1.7,C:-0.3 

# neb_make Usage: 
# Default: 
#     python neb_make.py [init_image] [final_image] [n_max] [optional]
# [optional]:
#     --fix [height]:[direction] : fix atom below height (fractional) in direction (0,1,2 for x,y,z)
#     --mag [element1]:[magmom1],[element2]:[magmom2],... : set initial magmom for atoms of element
# Use existing guess and do continuation calculation
#     python neb_make.py -i [init_guess_traj] [n_max] [optional]

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

# if plus dimer method
# sbatch dimer_submit.sh
