#!/bin/bash
#SBATCH -J Cu-Ag-NEB
#SBATCH -p amd41
##SBATCH -t 48:00:00
#SBATCH -n 16
#SBATCH -N 1
#SBATCH -o ase.log
#SBATCH -e ase.err

module load abacus/3.3.4-icc
module load miniconda3
conda activate ase

# OMP abacus have better performace than MPI - 2023.9.1
#export OMP_NUM_THREADS=4
#NP=`expr $SLURM_NTASKS / $OMP_NUM_THREADS`

#echo "change mincpus"
#scontrol update JobId=$SLURM_JOB_ID MinCPUsNode=4

# doing job
touch JobProcessing.state
echo `date` >> JobProcessing.state 

#mpirun -np $NP abacus
#python run_neb.py
export PYTHONPATH=$(realpath ../..):$PYTHONPATH
python run_neb.py

echo `date` >> $HOME/finish
echo `pwd` >> $HOME/finish
echo `date` >> JobProcessing.state
