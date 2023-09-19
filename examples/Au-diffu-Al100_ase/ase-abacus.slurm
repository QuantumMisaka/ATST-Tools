#!/bin/bash
#SBATCH -J ASE-NEB
#SBATCH -p amd41
##SBATCH -t 48:00:00
#SBATCH -n 64
#SBATCH -N 1
#SBATCH -o ase-abacus.log
#SBATCH -e ase-abacus.err

module load abacus
module load miniconda3
conda activate ase

#export OMP_NUM_THREADS=4
#NP=`expr $SLURM_NTASKS / $OMP_NUM_THREADS`

#echo "change mincpus"
#scontrol update JobId=$SLURM_JOB_ID MinCPUsNode=4

# doing job
touch JobProcessing.state
echo `date` >> JobProcessing.state 

python ./neb.py
#mpirun -np $NP abacus

echo `date` >> $HOME/finish
echo `pwd` >> $HOME/finish
echo `date` >> JobProcessing.state
