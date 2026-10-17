"""Real MPI regressions for API progress and peer-safe failure exits."""

from __future__ import annotations

import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[2]
MPI_TIMEOUT_SECONDS = 15


@pytest.mark.parametrize("mode", ["start_callback", "start_stream", "image_callback", "success"])
def test_two_rank_progress_reaches_a_bounded_common_outcome(tmp_path: Path, mode: str) -> None:
    """Root sink failures release peers, while successful events stay root-only."""
    if os.environ.get("ATST_RUN_MPI_TESTS") != "1":
        pytest.skip("set ATST_RUN_MPI_TESTS=1 to run real MPI launcher regressions")
    launcher = shutil.which("mpiexec") or str(Path(sys.executable).with_name("mpiexec"))
    if not Path(launcher).is_file():
        pytest.skip("mpiexec is unavailable")

    smoke = """
import io
import os
from mpi4py import MPI
from ase import Atoms
from atst_tools.api import RunOptions, run_workflow
from atst_tools.api import services
from atst_tools.api.models import WorkflowExecutionError

rank = MPI.COMM_WORLD.rank
mode = os.environ['ATST_PROGRESS_TEST_MODE']
events = []
class BrokenStream:
    def write(self, data):
        raise OSError('stream is unavailable')
    def flush(self):
        pass
def callback(event):
    events.append(event['event'])
    if mode == 'start_callback' and event['event'] == 'workflow_start':
        raise RuntimeError('start callback is unavailable')
    if mode == 'image_callback' and event['event'] == 'image_step':
        raise RuntimeError('image callback is unavailable')

dry_run = mode in ('start_callback', 'start_stream')
if not dry_run:
    services._dispatch_normalized = lambda *args: [Atoms('H')]
calculation = (
    {'type': 'relax', 'init_structure': 'unused.traj'} if dry_run else
    {'type': 'neb', 'init_chain': 'unused.traj'}
)
stream = BrokenStream() if mode == 'start_stream' else io.StringIO()
try:
    result = run_workflow(
        {'calculation': calculation,
         'calculator': {'name': 'abacus', 'abacus': {'parameters': {}}}},
        RunOptions(dry_run=dry_run, progress=True, progress_stream=stream,
                   progress_callback=callback),
    )
except WorkflowExecutionError as error:
    assert mode != 'success'
    if rank == 0:
        assert error.__cause__ is not None
        assert error.to_document()['cause']['type'] in ('RuntimeError', 'OSError')
    outcome = 1
else:
    assert mode == 'success'
    assert result.status == 'complete'
    assert events == (['workflow_start', 'image_step'] if rank == 0 else [])
    outcome = 0
assert MPI.COMM_WORLD.allreduce(outcome) == (0 if mode == 'success' else 2)
print('PROGRESS_MPI_OK', mode, rank, flush=True)
"""
    environment = dict(os.environ)
    environment["PYTHONPATH"] = str(ROOT / "src")
    environment["ATST_PROGRESS_TEST_MODE"] = mode
    process = subprocess.Popen(
        [launcher, "-n", "2", sys.executable, "-c", smoke],
        cwd=tmp_path,
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        start_new_session=True,
    )
    try:
        stdout, stderr = process.communicate(timeout=MPI_TIMEOUT_SECONDS)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        process.communicate()
        pytest.fail(f"MPI progress regression timed out after {MPI_TIMEOUT_SECONDS} seconds")
    assert process.returncode == 0, stdout + stderr
    assert stdout.count(f"PROGRESS_MPI_OK {mode}") == 2, stdout + stderr
