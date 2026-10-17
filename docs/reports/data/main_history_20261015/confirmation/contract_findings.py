"""Read-only CPU contract review reproductions for atst-tools HEAD 5ad5de5."""
import json
import os
from pathlib import Path
import runpy
import subprocess
import sys
import tempfile
import types
from unittest.mock import patch

ROOT = Path.cwd()
assert (ROOT / 'pyproject.toml').is_file(), 'Run this probe from the repository root'
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / 'src'))
from ase import Atoms
from ase.calculators.lj import LennardJones
from ase.io import write
from atst_tools.api import RunOptions, run_workflow
from atst_tools.calculators.factory import CalculatorFactory
import atst_tools.calculators.constant_potential as current
from atst_tools.utils.summary import summarize_artifact_manifest

print('HEAD:', subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip())
fixture = runpy.run_path(str(ROOT / 'tests/unit/test_constant_potential.py'))
old = types.ModuleType('cp_v228_review')
sys.modules[old.__name__] = old
exec(compile(subprocess.check_output(['git','-C',str(ROOT),'show','v2.2.8:src/atst_tools/calculators/constant_potential.py'],text=True), '<v2.2.8 constant_potential.py>', 'exec'), old.__dict__)
with tempfile.TemporaryDirectory(prefix='atst-contract-review-') as temp:
    a = Atoms('H', positions=[[0,0,0]], cell=[10]*3, pbc=True)
    a.calc = current.ConstantPotentialCalculator(fixture['_factory'], potential_v=1., reference_electrons=10., initial_electrons=8., work_ref=4.6, work_ref_source='test', capacitance_initial=1/.12, capacitance_unit='e/V', potential_tolerance_v=1e-12, directory=Path(temp)/'cp')
    a.get_potential_energy()
    payload = a.info[current.CP_FACTS_INFO_KEY]
    print('\nCP cross-version:')
    print('schema:',payload['schema'])
    print('current reader:',current._constant_potential_facts_error(a,payload))
    print('v2.2.8 reader:',old._constant_potential_facts_error(a,payload))
    print('current restart:',current.constant_potential_restart_initial_electrons(a,payload['identity'],1e-12))
    try:
        old.constant_potential_restart_initial_electrons(a,payload['identity'],1e-12)
    except old.ConstantPotentialError as error:
        print('v2.2.8 restart:',type(error).__name__,str(error))
    else:
        raise AssertionError('Expected v2.2.8 rejection')
    print('v2.2.8 allow_recompute:',old.constant_potential_restart_initial_electrons(a,payload['identity'],1e-12,allow_recompute=True))
    payload['facts']['profile_status'] = 'reference_only_unvalidated'
    print('legacy-field payload/current reader:',current._constant_potential_facts_error(a,payload))
    print('legacy-field payload/v2.2.8 reader:',old._constant_potential_facts_error(a,payload))

    print('\nNested manifest path:')
    original = Path.cwd()
    os.chdir(temp)
    try:
        chain = [Atoms('H2', positions=[[0,0,0],[distance,0,0]]) for distance in (1.05,1.10,1.15)]
        write('chain.traj', chain)
        config = {'calculator': {'name':'abacus','abacus':{'parameters':{}}}, 'calculation': {'type':'neb','init_chain':'chain.traj','parallel':False,'climb':False,'max_steps':1,'artifact_manifest':'reports/artifacts.json'}}
        with patch.object(CalculatorFactory, 'get_calculator', side_effect=lambda *args, **kwargs: LennardJones()):
            result = run_workflow(config)
        summary = summarize_artifact_manifest(result.artifact_manifest)
        print('workflow status:', result.status)
        print('manifest artifact:',result.artifacts)
        print('actual trajectory exists:',Path('neb.traj').exists())
        print('summary status:',summary['status'])
        print('summary artifacts:',summary['artifacts'])
        assert Path('neb.traj').exists()
        assert summary['status']['n_missing'] == 1
    finally:
        os.chdir(original)
