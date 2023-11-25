#!/bin/bash
#SBATCH -J FS
#SBATCH -p amd41
##SBATCH -t 48:00:00
#SBATCH -n 64
#SBATCH -N 1
#SBATCH -o abacus.log
#SBATCH -e abacus.err

module load abacus/3.4.1-icc

# OMP abacus have better performace than MPI - 2023.9.1
export OMP_NUM_THREADS=16
NP=`expr $SLURM_NTASKS / $OMP_NUM_THREADS`

#echo "change mincpus"
#scontrol update JobId=$SLURM_JOB_ID MinCPUsNode=4

# doing job
touch JobProcessing.state
echo `date` >> JobProcessing.state 

mpirun -np $NP abacus

echo `date` >> $HOME/finish
echo `pwd` >> $HOME/finish
echo `date` >> JobProcessing.state
