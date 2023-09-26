#!/bin/bash
export PYTHONPATH=$(realpath ../../../):$PYTHONPATH
mpirun -np 4 gpaw python run_neb.py | tee run_neb.log
