#!/bin/bash
#SBATCH -J AuPd-OH-CH3-H2O
#SBATCH -p C064M0256G
##SBATCH -t 48:00:00
#SBATCH --nodes=2
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=4
#SBATCH -o abacus.out
#SBATCH -e abacus.err

module load abacus/3.7.5-icx

export OMP_NUM_THREADS=4
#NP=`expr $SLURM_NTASKS / $OMP_NUM_THREADS`
NP=$SLURM_NTASKS

#echo "change mincpus"
#scontrol update JobId=$SLURM_JOB_ID MinCPUsNode=4

# doing job
touch JobProcessing.state
echo $SLURM_JOB_ID > JobProcessing.state
echo `date` >> JobProcessing.state 

srun hostname -s | sort -n > slurm.hosts
mpirun -machinefile slurm.hosts -np $NP abacus

echo `date` >> $HOME/finish
echo `pwd` >> $HOME/finish
echo `date` >> JobProcessing.state
