import numpy as np
import os, time, glob
from pathlib import Path
from typing import List
from dflow import(
    Step,
    Steps,
    Workflow,
    upload_artifact,
    download_artifact
)
from dflow.python import (
    OP,
    OPIO,
    Artifact,
    OPIOSign,
    PythonOPTemplate, 
    Parameter
)

from dflow.plugins.bohrium import BohriumContext, BohriumExecutor
# from dflow import config, s3_config
# config["host"] = "https://workflow.test.dp.tech"
# s3_config["endpoint"] = "39.106.93.187:30900"
# config["k8s_api_server"] = "https://60.205.59.4:6443"
# config["token"] = "eyJhbGciOiJSUzI1NiIsImtpZCI6IlhMRGZjbnNRemE4RGQyUXRMZG1MX3NXeG5TMzlQTnhnSHZkS1lGM25SODAifQ.eyJpc3MiOiJrdWJlcm5ldGVzL3NlcnZpY2VhY2NvdW50Iiwia3ViZXJuZXRlcy5pby9zZXJ2aWNlYWNjb3VudC9uYW1lc3BhY2UiOiJhcmdvIiwia3ViZXJuZXRlcy5pby9zZXJ2aWNlYWNjb3VudC9zZWNyZXQubmFtZSI6ImFyZ28tdG9rZW4tajd0a3MiLCJrdWJlcm5ldGVzLmlvL3NlcnZpY2VhY2NvdW50L3NlcnZpY2UtYWNjb3VudC5uYW1lIjoiYXJnbyIsImt1YmVybmV0ZXMuaW8vc2VydmljZWFjY291bnQvc2VydmljZS1hY2NvdW50LnVpZCI6IjBhNzI1N2JhLWZkZWQtNGI2OS05YWU2LTZhY2U0M2UxNjdlNiIsInN1YiI6InN5c3RlbTpzZXJ2aWNlYWNjb3VudDphcmdvOmFyZ28ifQ.Gg7pctEsZC-2ZkFjHv-q21mOzBeuThocTMoNV2ZaLtOuxvXQOiVQhS8nq8nyPBiygJ3okOXKPhhrEH8Oe0kWuYtEXc88e1kX_MarQLCXLYSN53cdTLlgZQn01hHaHLO6KJubgU8mymNKj260GjDSf35a7wt8NgQIwm9ftqEwYuPXrm2yZEnhtbuNgfdpLIhw_DQxLXvwjTiny7vwR7ANpHfaynf2l0E12il3C7xeTP-lcPUm9BSFObO3icUbz67n0qsz3j8QWxRdH-jTzIr7tTvFP8SpdJbvMBmI4fgU01FR5CnWx296I9bzXjbuNefZGNu9ZuJ5RLiQDt5xmbTweQ"
bohrium_context = BohriumContext(username="XXXXXX", password="XXXXXX",  executor="mixed",  extra={}     )
from dflow import config, s3_config
from dflow.plugins import bohrium
from dflow.plugins.bohrium import TiefblueClient
config["host"] = "https://workflows.deepmodeling.com"
config["k8s_api_server"] = "https://workflows.deepmodeling.com"
bohrium.config["username"] = "XXXXXX"
bohrium.config["password"] = "XXXXXX"
bohrium.config["project_id"] = "XXXXXX"
s3_config["repo_key"] = "oss-bohrium"
s3_config["storage_client"] = TiefblueClient()

def structure_optimization(structure_file: Path, pp: dict, basis: dict,trajectory_file: str) -> Path:
    r"""
    Optimize the structure using the Abacus calculator.

    Parameters
    ----------
    structure_file : `Path`
        Path to the input structure files.
    pp : `dict`
        Dictionary containing pseudopotential information for Abacus calculator.
    basis : `dict`
        Dictionary containing basis information for Abacus calculator.
    trajectory_file : `str`
        Path to the output trajectory file.

    Returns
    -------
    `Path`
        Path to the output trajectory file containing the optimized structure.
    """

    from ase.io import read, write
    from ase.constraints import FixAtoms
    from ase.optimize import BFGS
    from ase.calculators.abacus import Abacus, AbacusProfile
    from ase.neb import NEB
    from ase.neb import NEBTools
    import matplotlib.pyplot as plt

    # Load structure
    structure = read(structure_file,format="abacus")

    # Set up the calculator
    abacus='/usr/local/bin/abacus'
    profile = AbacusProfile(argv=['mpirun','-n','32',abacus])
    kpts = [2,2,2]
    structure.calc = Abacus(profile=profile, ntype=2, ecutwfc=100, scf_nmax=50, smearing_method='gaussian', smearing_sigma=0.01, basis_type='lcao',
                mixing_type='pulay', mixing_beta='0.7', scf_thr=1e-8, out_chg=1, calculation='scf', force_thr=0.001, stress_thr=5,
                cal_force=1, cal_stress=1, out_stru=1, pp=pp, basis=basis, kpts=kpts)

    # Optimize structure
    #mask = [atom.tag > 1 for atom in structure]
    #structure.set_constraint(FixAtoms(mask=mask))

    opt = BFGS(structure, trajectory=trajectory_file)
    opt.run(fmax=0.05)

    return Path(trajectory_file)

def neb_calculation(initial_neb_structure: Path, final_neb_structure: Path, num_images: int, pp: dict, basis: dict) -> float:
    r"""
    Perform the NEB (Nudged Elastic Band) calculation using the Abacus calculator.

    Parameters
    ----------
    initial_neb_structure : `Path`
        Path to the initial structure file for the NEB calculation.
    final_neb_structure : `Path`
        Path to the final structure file for the NEB calculation.
    num_images : `int`
        Number of interpolation images in the NEB calculation.
    pp : `dict`
        Dictionary containing pseudopotential information for Abacus calculator.
    basis : `dict`
        Dictionary containing basis information for Abacus calculator.

    Returns
    -------
    `float`
        The calculated barrier of the reaction.
    """

    from ase.io import read, write
    from ase.constraints import FixAtoms
    from ase.optimize import BFGS
    from ase.calculators.abacus import Abacus, AbacusProfile
    from ase.neb import NEB
    from ase.neb import NEBTools
    import matplotlib.pyplot as plt

    # Load initial and final structures
    initial = read(initial_neb_structure)
    final = read(final_neb_structure)

    # Set up the calculator
    abacus='/usr/local/bin/abacus'
    profile = AbacusProfile(argv=['mpirun','-n','32',abacus])
    kpts = [2,2,2]
    calc = Abacus(profile=profile, ntype=2, ecutwfc=100, scf_nmax=50, smearing_method='gaussian', smearing_sigma=0.01, basis_type='lcao',
                mixing_type='pulay', mixing_beta='0.7', scf_thr=1e-8, out_chg=1, calculation='scf', force_thr=0.001, stress_thr=5,
                cal_force=1, cal_stress=1, out_stru=1, pp=pp, basis=basis, kpts=kpts)

    # Make a band consisting of 'num_images' number of images
    constraint = FixAtoms(mask=[atom.tag > 1 for atom in initial])
    images = [initial]
    for i in range(num_images):
        image = initial.copy()
        image.calc = Abacus(profile=profile, ntype=2, ecutwfc=100, scf_nmax=50, smearing_method='gaussian', smearing_sigma=0.01, basis_type='lcao',
                mixing_type='pulay', mixing_beta='0.7', scf_thr=1e-8, out_chg=1, calculation='scf', force_thr=0.001, stress_thr=5,
                cal_force=1, cal_stress=1, out_stru=1, pp=pp, basis=basis, kpts=kpts)
        image.set_constraint(constraint)
        images.append(image)
    images.append(final)

    # Optimize NEB path
    neb = NEB(images)
    neb.interpolate('idpp')
    opt = BFGS(neb, trajectory='neb.traj')
    opt.run(fmax=0.05,steps=200)

    # Get the calculated barrier and the energy change of the reaction.
    images = read(f'neb.traj@-{num_images+2}:')
    nebtools = NEBTools(images)
    barrier, energy_change = nebtools.get_barrier()
    print(barrier)
    return barrier

class StructureOptimization(OP):
    @classmethod
    def get_input_sign(cls) -> OPIOSign:
        return OPIOSign({
            "opt_input_files": Artifact(Path),
            "pp": Parameter(dict),
            "basis": Parameter(dict)
        })

    @classmethod
    def get_output_sign(cls) -> OPIOSign:
        return OPIOSign({
            "opt_output_files": Artifact(List[Path]),
        })

    @OP.exec_sign_check
    def execute(self, op_in: OPIO) -> OPIO:
        os.system("pip install git+https://gitlab.com/1041176461/ase-abacus.git")
        initial_structure = op_in["opt_input_files"]/"initial_stru"
        final_structure = op_in["opt_input_files"]/"final_stru"
        pp = op_in["pp"]
        basis = op_in["basis"]

        os.chdir(op_in["opt_input_files"])
        initial_opt_traj = structure_optimization(initial_structure,pp,basis,"initial.traj")
        final_opt_traj = structure_optimization(final_structure,pp,basis,"final.traj")

        op_out = {
            "opt_output_files": [initial_opt_traj, 
                                 final_opt_traj, 
                                 op_in["opt_input_files"]/"Li_ONCV_PBE-1.2.upf", 
                                 op_in["opt_input_files"]/"Si_ONCV_PBE-1.2.upf",
                                 op_in["opt_input_files"]/"Li_gga_8au_100Ry_4s1p.orb",
                                 op_in["opt_input_files"]/"Si_gga_8au_100Ry_2s2p1d.orb"]
        }
        return op_out

class NEBCalculation(OP):
    @classmethod
    def get_input_sign(cls) -> OPIOSign:
        return OPIOSign({
            "neb_input_files": Artifact(Path),
            "num_images": Parameter(int),
            "pp": Parameter(dict),
            "basis": Parameter(dict)
        })

    @classmethod
    def get_output_sign(cls) -> OPIOSign:
        return OPIOSign({
            "barrier": Parameter(float)
        })

    @OP.exec_sign_check
    def execute(self, op_in: OPIO) -> OPIO:
        os.system("pip install git+https://gitlab.com/1041176461/ase-abacus.git")

        initial_neb_structure = op_in["neb_input_files"]/"initial.traj"
        final_neb_structure = op_in["neb_input_files"]/"final.traj"

        num_images = op_in["num_images"]
        pp = op_in["pp"]
        basis = op_in["basis"]

        os.chdir(op_in["neb_input_files"])
        barrier = neb_calculation(initial_neb_structure,final_neb_structure, num_images, pp, basis=basis)

        op_out = {
            "barrier": barrier
        }
        return op_out


if __name__ == "__main__":

    artifact0 = upload_artifact(["./initial_stru",
                                 "./final_stru",
                                 "./Li_ONCV_PBE-1.2.upf",
                                 "./Si_ONCV_PBE-1.2.upf",
                                 "./Li_gga_8au_100Ry_4s1p.orb",
                                 "./Si_gga_8au_100Ry_2s2p1d.orb"])

    pp = {'Li':'Li_ONCV_PBE-1.2.upf','Si':'Si_ONCV_PBE-1.2.upf'}
    basis = {'Li':'Li_gga_8au_100Ry_4s1p.orb','Si':'Si_gga_8au_100Ry_2s2p1d.orb'}

    # Create steps
    structure_optimization_step = Step(
        name="structure-optimization",
        template=PythonOPTemplate(StructureOptimization,image="registry.dp.tech/dptech/abacus:3.2.1"),
        executor=BohriumExecutor(executor="bohrium_v2", extra={"scassType":"c64_m64_cpu","projectId":11176, "jobType":"container"}),
        artifacts={"opt_input_files":artifact0},
        parameters={"pp": pp, "basis": basis}
    )

    neb_calculation_step = Step(
        name="neb-calculation",
        template=PythonOPTemplate(NEBCalculation,image="registry.dp.tech/dptech/abacus:3.2.1"),
        executor=BohriumExecutor(executor="bohrium_v2", extra={"scassType":"c64_m64_cpu","projectId":11176, "jobType":"container"}),
        artifacts={"neb_input_files": structure_optimization_step.outputs.artifacts["opt_output_files"]},
        parameters={"num_images": 7,"pp": pp, "basis": basis}
    )

    # Create workflow
    wf = Workflow(name="neb-workflow", context=bohrium_context, host="https://workflow.test.dp.tech")
    wf.add(structure_optimization_step)
    wf.add(neb_calculation_step)

    # Submit workflow
    wf.submit()

