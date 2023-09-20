# ase-neb-abacus
NEB or CI-NEB method by ASE and ABACUS using ASE-ABACUS interfaces

## Dependencies:
- ASE
- ASE-ABACUS interfaces
- ABACUS

## Usage
- `abacus_neb.py` have all the NEB workflow by using ASE-ABACUS interface to do NEB calculation
- examples can be run by `python run_neb.py`, or simply bt `./run.sh`
- Before use, make `abacus_neb.py` detectable by python, e.g. `export PYTHONPATH=/path/to/ase-neb-abacus:$PYTHONPATH`

## Method
- For serial NEB calculation, DyNEB, namely dynamic NEB method is for default used.
- Users can change lots of parameter for different NEB setting. one can refer to [ASE NEB calculator](https://wiki.fysik.dtu.dk/ase/ase/neb.html#module-ase.neb) for more details: 

## Examples
- Au-diffu-Al100_ase: Au diffusion on Al(100) surface, an example for running ASE-NEB-ABACUS based on ASE-constructed initial and final states
- Ag-diffu-Cu100: Ag diffusion on Cu(100) surface, an example for running ASE-NEB-ABACUS based on ABACUS-calculated initial and final states. 
- Li-diffu-Si: Li diffusion in Si, an example for running ASE-NEB-ABACUS based on existing ABACUS input files of initial and final state, using ABACUS as SCF calculator and ASE as optimizer and NEB calculator.  Also, an dflow example is proposed.
- N2-Cu111: N2 dissociation on Cu(111) surface, an example for running ASE-NEB-ABACUS based on existing ABACUS input files of initial and final state, use ABACUS as optimizer for initial and final state and use ASE as NEB calculator. 

## Developing
- [x] Use interface to read ABACUS STRU file and ABACUS output
- [x] Flexible input for different NEB method in ASE
- [ ] More test in surface reaction system
- [ ] Parallel computing during images relaxation


