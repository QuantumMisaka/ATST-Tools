#!/bin/bash
#SBATCH -J Pt111-CO
#SBATCH -p C064M0256G
#SBATCH -t 48:00:00
#SBATCH -n 64
#SBATCH -N 1
#SBATCH -o abacus.log
#SBATCH -e abacus.err

module load abacus/3.4.1-icx

export OMP_NUM_THREADS=4
NP=`expr $SLURM_NTASKS / $OMP_NUM_THREADS`

#echo "change mincpus"
#scontrol update JobId=$SLURM_JOB_ID MinCPUsNode=4

# doing job
touch JobProcessing.state
echo $SLURM_JOB_ID >> JobProcessing.state
echo `date` >> JobProcessing.state 

mpirun -np $NP abacus

echo `date` >> $HOME/finish
echo `pwd` >> $HOME/finish
echo `date` >> JobProcessing.state
