# ATST-Tools
Advanced ASE Transition State Tools for ABACUS, including:
- NEB and improvement, like Dynamic NEB
- CI-NEB, IT-NEB and others
- AutoNEB
- Dimer

Version v1.2.1

## Dependencies:
- [ASE](https://wiki.fysik.dtu.dk/ase/about.html)
- [ABACUS](https://abacus.deepmodeling.com/en/latest/)
- [ASE-ABACUS interface](https://gitlab.com/1041176461/ase-abacus)
- [GPAW](https://wiki.fysik.dtu.dk/gpaw/install.html) if one wants to run NEB images relaxation in parallel

Notice: GPAW and ABACUS should be dependent on same MPI and libraries environments. 
For instance, if your ABACUS is installed by Intel-OneAPI toolchain, your GPAW should NOT be dependent on gcc-toolchain like OpenMPI and OpenBLAS.

ATST-Tools is Under actively development, please let me know if any problem occurs.

## Workflow

![ATST-workflow](./img/ATST-workflow.png)


## Developing
- [x] Use interface to read ABACUS STRU file and ABACUS output
- [x] Flexible input for different NEB method in ASE
- [x] `DyNEB` implementation and test
- [x] Now used optimum option: idpp-guess + IT-NEB + CI-NEB parallel method
- [x] Give an initial guess print-out of NEB images
- [x] Decoupling to init-guess -> NEB calculation -> result post-process
- [x] Parallel computing during images relaxation by `gpaw python`
- [x] `AutoNEB` implementation
- [x] Give user a way to specify constraint information (along one direction) and initial magnetic moments (for each atom of one element) for NEB initlal guess chain and two end points (while these information cannot be read from calculation result)
- [x] More test in surface reaction system
- [x] Connected to Dimer method
- [x] More test in magnetic surface reaction system
- [x] Put calculation setting in an independent file (decoupling *_run.py)
- [ ] More fiexible options for NEB, Dimer and AutoNEB, like full properties in trajectory file, and fiexibly utilize SCF wavefunction/charge output files from previous calculation.

## Workflow libraries files and setting
All workflow libraries files and re-constructed ASE libraries will be put in `./source` directory. including:
```bash
source
├── abacus_autoneb.py
├── abacus_dimer.py
└── abacus_neb.py
```
Before use running scripts, you should add these libraries into your PYTHONPATH:
```bash
export PYTHONPATH=/path/to/source:$PYTHONPATH
```

## NEB workflow

### Method
- For serial NEB calculation, DyNEB, namely dynamic NEB method `ase.mep.neb.DyNEB` is for default used.
- For parallel NEB calculation, `ase.mep.neb.NEB` traditional method is for default used.
- The Improved Tangent NEB method `IT-NEB` and Climbing Image NEB method `CI-NEB` in ASE are also default used in this workflow, which is highly recommended by Sobervea. In `AutoNEB`, `eb` method is used for default, but Improved Tangent method is also recommended.
- Users can change lots of parameter for different NEB setting. one can refer to [ASE NEB calculator](https://wiki.fysik.dtu.dk/ase/ase/neb.html#module-ase.neb) for more details: 
- The workflow also support use `AutoNEB` method in ASE。You can view AutoNEB method in paper below. Also. One can refer to [AutoNEB](https://wiki.fysik.dtu.dk/ase/ase/neb.html#autoneb) to view it. 
> E. L. Kolsbjerg, M. N. Groves, and B. Hammer, J. Chem. Phys, 145, 094107, 2016. (doi: 10.1063/1.4961868)

The AutoNEB method in ASE lies in `ase.mep.autoneb.AutoNEB` object, which will do NEB calculation in following steps:
1. Define a set of images and name them sequentially. Must at least have a relaxed starting and ending image. User can supply intermediate guesses which do not need to have previously determined energies (probably from another
NEB calculation with a lower level of theory)
1. AutoNEB will first evaluate the user provided intermediate images
2. AutoNEB will then add additional images dynamically until n_max is reached
3. A climbing image will attempt to locate the saddle point
4. All the images between the highest point and the starting point are further relaxed to smooth the path
5. All the images between the highest point and the ending point are further relaxed to smooth the path

Step 4 and 5-6 are optional steps. Note that one can specify different `fmax` for CI-NEB in step 4 compared with other NEB calculation, which can be set as `fmax=[fmax1, fmax2]` in `AutoNEB` object.

> Notice: in surface calculation and hexangonal system, the vaccum and c-axis should be set along y-direction but not z-direction, which is much more efficient for ABACUS calculation.


### Usage
#### Basic NEB
The NEB workflow is based on 3 main python scripts and 1 workflow submit script. Namely:

- `neb_make.py` will make initial guess for NEB calculation, which is based on ABACUS (and other calculator) output files of initial and final state. This script will generate `init_neb_chain.traj` for neb calculation. Also, You can do continuation calculation by using this script. You can get more usage by `python neb_make.py`. 
- `neb_run.py` is the key running script of NEB, which will run NEB calculation based on `init_neb_chain.traj` generated by `neb_make.py`. This script will generate `neb.traj` for neb calculation. Users should edit this file to set parameters for NEB calculation. sereal running can be performed by `python neb_run.py`, while parallel running can be performed by `mpirun gpaw python neb_run.py`.
When running, the NEB trajectory will be output to `neb.traj`, and NEB images calculation will be doing in `NEB-rank{i}` directory for each rank which do calculation of each image. 
- `neb_post.py` will post-process the NEB calculation result, which will based on `neb.traj` from neb calculation. This script will generate nebplots.pdf to view neb calculation result, and also print out the energy barrier and reaction energy. You can get more usage by `python neb_post.py`. Meanwhile, users can also view result by `ase -T gui neb.traj` or `ase -T gui neb.traj@-{n_images}:` by using ASE-GUI
- `neb_submit.sh` will do all NEB process in one workflow scripts and running NEB calculation in parallel. Users should edit this file to set parameters for NEB calculation. Also this submit script can be used as a template for job submission in HPC. the Default setting is for `slurm` job submission system.

#### AutoNEB method 
In ATST-Tools, the AutoNEB method can be easily used by the following scripts
- `autoneb_run.py` is the key running script for `AutoNEB` method, which is like `neb_run.py` but the NEB workflow in `AutoNEB` is enhanced and the I/O logic have some difference. Users can use it with `mpirun gpaw python autoneb_run.py` by existing `init_neb_chain.traj` which can only contain initial and final state or contain some initial-guess.
- `autoneb_submit.sh` will do all AutoNEB process in one workflow and running AutoNEB calculation in parallel. Users should edit this file to set parameters for AutoNEB calculation. Also this submit script can be used as a template for job submission in HPC. the Default setting is for `slurm` job submission system.
- `neb_make.py` and `neb_post.py` can be used for `AutoNEB` method, but the workflow have slight difference. 

#### Running
Users can run NEB each step respectively: 
1. `python neb_make.py [INIT/result] [FINAL/result] [n_max]` to create initial guess of neb chain
   1. Also You can use `python neb_make.py -i [input_traj_file] [n_max]` to create initial guess from existing traj file, which can be used for continuation calculation.
2. `python neb_run.py` or `mpirun -np [nprocs] gpaw python neb_run.py` to run NEB calculation
3. `python neb_post.py neb.traj [n_max]` to post-process NEB calculation result

Users can run AutoNEB each step respectively:
1. `python neb_make.py [INIT/result] [FINAL/result] [nprocs]` to create initial guess of neb chain
2. `mpirun -np [nprocs] gpaw python autoneb_run.py` to run AutoNEB calculation
3. `python neb_post.py --autoneb run_autoneb???.traj` to post-process NEB calculation result

Also, user can run each step in one script `neb_submit.sh` by `bash neb_submit.sh` or `sbatch neb_submit.sh`. AutoNEB scripts usage is like that. 

> Notice: Before you start neb calculation process, make sure that you have check the nodes and cpus setting and other setting like n_max, constraints and initial magnetic moments in `*neb_submit.sh` and `*neb_run.py` to make sure that you can reach the highest performance and reach the simulation result you want !!!   


#### Continuation calculation for NEB
If NEB or AutoNEB is break down somehow, you can do continuation calculation based on saved trajectory files and ATST-Tools scripts.

For NEB, you can simply:
```bash
python neb_make.py -i neb.traj [n_max]
```
to generate `init_neb_chain.traj` for continuation calculation. You can also `python neb_post.py neb.traj` to generate the latest neb band `neb_latest.traj` and do continuation calculation by `python neb_make.py -i neb_latest.traj [n_max]`. note that `n_max = n_image - 2`

For AutoNEB, you need to get `neb_latest.traj` in a more compicated way:
```bash
python neb_post.py --autoneb ./AutoNEB_iter/run_autoneb???iter[index].traj
```
to generate `neb_latest.traj` from certain index (like 006) stage of AutoNEB calculation, and then
```bash
python neb_make.py -i neb_latest.traj [n_max]
```
to generate `init_neb_chain.traj` for continuation calculation.

> Note: Linux shell will automatically detect and sort number of index, so you will not be worried about using format like `run_autoneb???iter005.traj`, the consequence will be right, for example:
```bash
 test> ll run_autone*.traj
-rw-r--r-- 1 james james 6.0K Nov 24 20:35 run_autoneb000.traj
-rw-r--r-- 1 james james 6.4K Nov 24 20:35 run_autoneb001.traj
-rw-r--r-- 1 james james 6.4K Nov 24 20:35 run_autoneb002.traj
-rw-r--r-- 1 james james 531K Nov 24 20:35 run_autoneb003.traj
-rw-r--r-- 1 james james 531K Nov 24 20:35 run_autoneb004.traj
-rw-r--r-- 1 james james 531K Nov 24 20:35 run_autoneb005.traj
-rw-r--r-- 1 james james 531K Nov 24 20:35 run_autoneb006.traj
-rw-r--r-- 1 james james 6.4K Nov 24 20:35 run_autoneb007.traj
-rw-r--r-- 1 james james 6.4K Nov 24 20:35 run_autoneb008.traj
-rw-r--r-- 1 james james 6.5K Nov 24 20:35 run_autoneb009.traj
-rw-r--r-- 1 james james 6.5K Nov 24 20:35 run_autoneb010.traj
-rw-r--r-- 1 james james 6.5K Nov 24 20:35 run_autoneb011.traj
-rw-r--r-- 1 james james 6.5K Nov 24 20:35 run_autoneb012.traj
-rw-r--r-- 1 james james 531K Nov 24 20:37 run_autoneb025.traj
```

#### Other scripts
Because ATST is originally based on ASE, the trajectory file can be directly read, view and analysis by `ase gui` and other ASE tools. Abide by `neb_make.py` and `neb_post.py`, We also offer some scripts to help you:
- `neb_dist.py`: This script will give distance between initial and final state, which is good for you to check whether the atoms in two image is correspondent, and is also a reference for setting number of n_max
- `traj_transform.py`: This script can transfer traj files into other format like `extxyz`, `abacus`(STRU), `cif` and so on (coming soon). Also if user specify `--neb` option, this script will automatically detect and cut the NEB trajectory when doing format transform. This script will be helpful for analysis and visualization of NEB trajectory.


## Dimer workflow

### Method
(Waiting for update)


### Usage
The Dimer workflow is based on 2 main python scripts and 2 workflow submit script. Namely:
- `neb2dimer.py` can be used by `python neb2dimer [neb.traj] ([n_max])`, which will transform NEB trajetory `neb.traj` or NEB result trajectory `neb_result.traj` to Dimer input files,  including:
- - `dimer_init.traj` for initial state of Dimer calculation, which is the highest energy image, namely, TS state. 
- - `displacement_vector.npy` for displacement vector of Dimer calculation, which will be generated from position minus of the nearest image before and after TS point, and be normalized to 0.01 Angstrom. 
- `dimer_run.py` is the key running script of Dimer calculation, which will run Dimer calculation based on `dimer_init.traj` and `displacement_vector.npy` generated by `neb2dimer.py` or based on other setting. This script will generate `dimer.traj` for Dimer calculation trajectory. Users should edit this file to set parameters for Dimer calculation, and run Dimer calculation by `python dimer_run.py`. When running, any Dimer images calculation will be doing in `Dimer` directory.
- `dimer_submit.sh` will do Dimer workflow in one scripts. The Default setting is for `slurm` job submission system.
- `neb-dimer_srun.sh` is a try to run NEB + Dimer calculation in one server scripts. The Default setting is for `slurm` job submission system.


## Notices
### Property Loss in Trajectory 
Some property should be get via specific way from trajectory files, and some will be lost in trajetory files, 
- Stress property will not be stored in trajetory file
- In NEB calculation, the Force property for fixed atoms and Stress property will NOT be stored in trajectroy file.
- in Dimer calculation, the Energy, Forces and Stress property will NOT be stored in trajetory file.
- in AutoNEB calculation, all property in processing trajectory will be stored in AutoNEB_iter directory, but in the result `run_autoneb???.traj`, the forces and stress information will be lost.


## Examples
- Li-diffu-Si: Li diffusion in Si, very easy example for serial and parallel NEB calculation
- H2-Au111: H2 dissociation on Au(111) surface. which will have NEB, AutoNEB and Dimer example. The barrier is around 1.1 eV consistent with existing paper and calculation result.
- N2-Cu111 : N2 dissociation on Cu(111) surface. which have high barrier as 4.0 eV but it's not hard to calculate
- CO-Pt111 : CO dissociation on Pt(111) surface. which have high barrier as 3.6 eV and including diffusion part. AutoNEB will always fail due to the low starting image. 
- Cy-Pt_graphene: Cyclohexane dehydrogenation on Pt-doped graphene surface. The barrier is around 1.3 eV. Noted that the `IT-NEB` result is wrong, but which is consistent to the result in VTST-Tools when using 4 image to do IT-NEB calculation.
- C-C coupling on Fe5C2(510) surface: C-C coupling on Fe5C2(510) surface with CO co-adsorption, which have 0.48 eV barrier in calculation result and existing paper. This is a typical example for reaction on magnetic surfaces.

More examples is welcomed from users. 