"""Assert reviewed defects using CPU fixtures without changing repository files.

Each invocation creates a fresh directory under /tmp for every workflow artifact.
EMT replaces the external calculator; the Sella IRC completion fixture isolates
post-processing, and an IPOPT shell permits real PyDMF path construction only.
No ABACUS, DeePMD, Slurm, or IPOPT solve is executed. Assertions describe current
bugs: a repaired implementation should make this evidence probe fail.
"""
from pathlib import Path
import os
import sys
import tempfile
import types
import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.calculators.emt import EMT
from atst_tools.workflows import execution, irc
from atst_tools.utils.restart_helpers import get_last_neb_band

artifact_dir = Path(tempfile.mkdtemp(prefix='atst-science-findings-', dir='/tmp'))
os.chdir(artifact_dir)
print('ARTIFACT_DIRECTORY', artifact_dir)
execution.CalculatorFactory.get_calculator = staticmethod(lambda *a, **kw: EMT())


def expect_exception(operation, expected_type, message_fragment, label):
    """Assert that the current defect raises the expected exception."""
    try:
        operation()
    except expected_type as exc:
        assert message_fragment in str(exc), (label, str(exc))
        print(label, type(exc).__name__, str(exc))
    else:
        raise AssertionError(f'{label}: expected {expected_type.__name__}')


# Actual initial NEB and restart: only the external calculator is replaced.
for name, x in [('start', 0), ('ts', .3), ('end', .6)]:
    write(name + '.traj', Atoms('H', positions=[[x, 0, 0]]))
neb_config = {'type': 'neb', 'trajectory': 'actual-neb.traj', 'parallel': False,
    'fmax': 1., 'max_steps': 1, 'make': {'n_images': 3,
        'init_structure': 'start.traj', 'final_structure': 'end.traj',
        'ts_guess': 'ts.traj', 'method': 'linear', 'no_align': True}}
actual_band = execution.run_neb({}, 'abacus', neb_config)
assert len(actual_band) == 6
recorded_band = read('actual-neb.traj', index=':')
assert len(recorded_band) == 6
print('PUBLIC NEB initial run completed, actual band size:', len(actual_band),
      'restart expected band size:', neb_config['make']['n_images'] + 2)
neb_config['restart'] = True
expect_exception(lambda: execution.run_neb({}, 'abacus', neb_config), ValueError,
    'not a whole number of bands with 5 images', 'PUBLIC NEB restart exception:')

# Five valid six-frame bands make the incorrect five-frame divisibility pass.
write('five-bands.traj', recorded_band * 5)
assert len(read('five-bands.traj', index=':')) == 30
selected = get_last_neb_band('five-bands.traj', 5)
assert len(selected) == 5
assert np.isclose(selected[0].positions[0, 0], .15)
neb_config['trajectory'] = 'five-bands.traj'
restarted_band = execution.run_neb({}, 'abacus', neb_config)
assert len(restarted_band) == 5
assert np.isclose(restarted_band[0].positions[0, 0], .15)
assert np.isclose(recorded_band[0].positions[0, 0], 0.)
assert not np.allclose(restarted_band[0].positions, recorded_band[0].positions)
print('PUBLIC NEB 30-frame restart completed, actual initial x:',
      restarted_band[0].positions[0, 0], 'expected initial x:',
      recorded_band[0].positions[0, 0], 'actual band size:', len(restarted_band))

# Real ASE FIRE descent writes result-bearing branch files but strips main files.
np.save('mode.npy', np.array([[1., 0, 0]]))
irc.IRCWorkflow({}, 'abacus', {'type': 'irc', 'backend': 'descent',
    'init_structure': 'start.traj', 'mode_vector': 'mode.npy',
    'max_steps': 1, 'trajectory': 'irc.traj',
    'normalized_trajectory': 'norm_irc.traj', 'direction': 'both'}).run()
for filename in ['irc_forward.traj', 'irc_reverse.traj', 'irc.traj', 'norm_irc.traj']:
    frames = read(filename, index=':')
    has_results = [a.calc is not None for a in frames]
    if filename in ['irc_forward.traj', 'irc_reverse.traj']:
        assert len(frames) == 1 and all(has_results)
        assert np.isfinite(frames[0].get_potential_energy())
        assert np.isfinite(frames[0].get_forces()).all()
    else:
        assert len(frames) == 2 and not any(has_results)
    print('IRC descent output:', filename, 'calculators:', has_results)

# A completed IRC optimizer isolates real workflow schema, IO, and normalization.
import sella
completed_directions = []
class CompletedIRC:
    def __init__(self, atoms, trajectory, **kwargs):
        self.atoms = atoms
        self.trajectory = trajectory
        self.nsteps = 0
    def run(self, fmax, steps, direction):
        self.atoms.get_potential_energy()
        self.atoms.get_forces()
        self.trajectory.write(self.atoms)
        completed_directions.append(direction)
        print('IRC controlled optimizer completed:', direction)
        return True
sella.IRC = CompletedIRC
workflow = irc.IRCWorkflow({}, 'abacus', {'type': 'irc',
    'init_structure': 'start.traj', 'trajectory': 'completed-sella-irc.traj'})
assert workflow.direction == 'both' and workflow.normalized_traj_file is None
expect_exception(workflow.run, TypeError, 'NoneType', 'PUBLIC IRC post-run exception:')
assert completed_directions == ['forward', 'reverse']

# Stub only IPOPT's shell; use the real DMF workflow validation/path constructor.
class Problem:
    def __init__(self, *a, **kw): pass
    def add_option(self, *a, **kw): pass
sys.modules['cyipopt'] = types.SimpleNamespace(Problem=Problem)
from atst_tools.workflows.dmf import DMFWorkflow
from atst_tools.external.pydmf.dmf import DirectMaxFlux
start = Atoms('HO', positions=[[0, 0, 0], [1, 0, 0]])
end = Atoms('OH', positions=[[0, 0, 0], [2, 0, 0]])
workflow = DMFWorkflow({}, 'abacus', {'type': 'dmf', 'init_file': 'i',
    'final_file': 'f', 'initial_path': 'linear',
    'remove_rotation_and_translation': False})
workflow._get_calc = lambda index: None
workflow._validate_pbc(start, end)
dmf = workflow._build_initial_path([start, end], DirectMaxFlux, None, {})
assert np.array_equal(dmf.images[-1].numbers, start.numbers)
assert not np.array_equal(dmf.images[-1].numbers, end.numbers)
np.testing.assert_allclose(dmf.images[-1].positions, end.positions)
print('DMF requested endpoint:', end.get_chemical_symbols(), end.positions.tolist())
print('DMF actual endpoint:', dmf.images[-1].get_chemical_symbols(),
      dmf.images[-1].positions.tolist())
print('ALL FOUR CURRENT DEFECTS CONFIRMED BY ASSERTIONS')
