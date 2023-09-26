# ase-neb-abacus
NEB or CI-NEB workflow by ASE and ABACUS using ASE-ABACUS interfaces,
Which use enhanced NEB method in ASE, like Dynamic NEB `DyNEB` or `AutoNEB`.

Version v0.3.1

## Dependencies:
- [ASE](https://wiki.fysik.dtu.dk/ase/about.html)
- [ABACUS](https://abacus.deepmodeling.com/en/latest/)
- [ASE-ABACUS interfaces](https://gitlab.com/1041176461/ase-abacus)
- [GPAW](https://wiki.fysik.dtu.dk/gpaw/install.html) if one wants to run NEB images relaxation in parallel

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
- N2-Cu111: N2 dissociation on Cu(111) surface, an example for running ASE-NEB-ABACUS based on existing ABACUS input files of initial and final state, use ABACUS as SCF calculator, use ASE for initial and final state optimization and use ASE as NEB calculator. 
- H2-Au111: H2 dissociation on Au(111) surface, an example for running ASE-NEB-ABACUS based on existing ABACUS input and output files of initial and final state, use ABACUS as optimizer for optimization of initial state and final state, use ASE as NEB calculator.
- CH4-Ir001: CH4 dissociation on Ir(001) surface, an example for running ASE-NEB-ABACUS based on existing ABACUS input files of initial and final state, use ABACUS as SCF calculator, use ASE for initial and final state optimization and use ASE as NEB calculator.
- CO-Pt111: CO dissociation on Pt(111) surface, an example for running ASE-NEB-ABACUS based on existing ABACUS input files of initial and final state, use ABACUS as SCF calculator, use ASE for initial and final state optimization and use ASE as NEB calculator. 
- Cy-Pt_graphene: Cyclohexane dehydrogenation on Pt-doped graphene surface, an example for running ASE-NEB-ABACUS based on existing ABACUS input files of initial and final state, use ABACUS as optimizer for optimization of initial state and final state, use ASE as NEB calculator


## Next Examples
- K-diffu-Fe5C2-510 : K diffusion on Fe5C2(510) surface, ,in this example we will focus on setting `magmom` during NEB calculation, which is important for spin-polarized magnetic system.
- K2O-diffu-Fe5C2-510: likely as K-diffu-Fe5C2-510
- CO-Fe100: CO dissociation on Fe(100) surface
- CO-Fe5C2-510: CO dissociation on Fe5C2(510) surface. Which is the final goal.
- `AutoNEB` implementation and test 

## Developing
- [x] Use interface to read ABACUS STRU file and ABACUS output
- [x] Flexible input for different NEB method in ASE
- [x] `DyNEB` implementation and test
- [x] Now used optimum option: idpp + DyNEB + IT-NEB + CI-NEB method
- [x] Make bottom atom fixed when read from `running*.log` of ABACUS
- [x] Give an initial guess print-out of NEB images
- [x] More test in surface reaction system
- [ ] Parallel computing during images relaxation by `gpaw python`
- [ ] More test in magnetic surface reaction system
- [ ] `AutoNEB` implementation


