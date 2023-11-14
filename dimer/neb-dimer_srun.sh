#!/bin/bash

# full workflow controling NEB and Dimer calculation

NEB_JOBID=$(sbatch neb_submit.sh) | awk '{print $4}' # get the jobid of neb and submit neb job
echo "NEB_JOBID is ${NEB_JOBID}"
DIMER_JOBID=$(sbatch --dependency=afterok:${NEB_JOBID} dimer_submit.sh) | awk '{print $4}' # submit dimer job after neb job done
echo "Dimer Job Submitted and will be run after NEB job !"
echo "DIMER_JOBID is ${DIMER_JOBID}"