#!/bin/bash
export PYTHONPATH=$(realpath ../../../):$PYTHONPATH
python run_neb.py | tee run_neb.log
