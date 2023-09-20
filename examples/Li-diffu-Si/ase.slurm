#!/bin/bash
#SBATCH -J Li-Si_diffu
#SBATCH -p C064M0256G
#SBATCH -t 48:00:00
#SBATCH -n 64
#SBATCH -N 1
#SBATCH -o ase.log
#SBATCH -e ase.err

module load abacus/3.3.4-icc
source /lustre/home/2201110432/apps/miniconda3/etc/profile.d/conda.sh
conda activate ase

#export OMP_NUM_THREADS=16
#NP=`expr $SLURM_NTASKS / $OMP_NUM_THREADS`

#echo "change mincpus"
#scontrol update JobId=$SLURM_JOB_ID MinCPUsNode=4

# doing job
touch JobProcessing.state
echo $SLURM_JOB_ID >> JobProcessing.state
echo `date` >> JobProcessing.state 

export PYTHONPATH=$(realpath ..):$PYTHONPATH
python run_neb.py
#mpirun -np $NP abacus

echo `date` >> $HOME/finish
echo `pwd` >> $HOME/finish
echo `date` >> JobProcessing.state
