"""Paired CLI/API checks for the shared scientific execution layer."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.emt import EMT
from ase.io import read, write

from helpers import FakeWorld, child_test_environment


def _write_relax_fixture(directory: Path) -> Path:
    """Write identical YAML and a periodic ASE structure for one entry point."""
    directory.mkdir()
    atoms = Atoms(
        "H2",
        positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.7]],
        cell=[10.0, 11.0, 12.0],
        pbc=[True, False, True],
    )
    write(directory / "initial.traj", atoms)
    config = directory / "config.yaml"
    config.write_text(
        "calculation:\n"
        "  type: relax\n"
        "  init_structure: initial.traj\n"
        "  trajectory: relax.traj\n"
        "  logfile: relax.log\n"
        "  fmax: 10.0\n"
        "  max_steps: 1\n"
        "calculator:\n"
        "  name: abacus\n"
        "  abacus:\n"
        "    parameters: {}\n",
        encoding="utf-8",
    )
    return config


def test_cli_and_api_execute_same_relax_fixture(monkeypatch, tmp_path):
    """Both entry points preserve scientific results and their manifest policy."""
    from atst_tools.api import RunOptions, run_workflow
    from atst_tools.scripts import main as cli_main
    from atst_tools.workflows import relax

    monkeypatch.setattr(
        relax.CalculatorFactory,
        "get_calculator",
        lambda *args, **kwargs: EMT(),
    )
    cli_dir = tmp_path / "cli"
    api_dir = tmp_path / "api"
    _write_relax_fixture(cli_dir)
    _write_relax_fixture(api_dir)

    monkeypatch.chdir(cli_dir)
    cli_result = cli_main.main(["config.yaml"])
    cli_atoms = read("final_relaxed.traj")
    cli_manifest = json.loads(Path("atst_artifacts.json").read_text())

    monkeypatch.chdir(api_dir)
    api_result = run_workflow("config.yaml", RunOptions(world=FakeWorld()))
    api_atoms = read("final_relaxed.traj")
    api_manifest = json.loads(Path("atst_artifacts.json").read_text())

    assert cli_result.status == api_result.status == "complete"
    assert cli_result.final_atoms is None
    assert api_result.final_atoms is not None
    for actual in (cli_atoms, api_atoms, api_result.final_atoms):
        np.testing.assert_allclose(actual.cell.lengths(), [10.0, 11.0, 12.0])
        assert actual.pbc.tolist() == [True, False, True]
    np.testing.assert_allclose(cli_atoms.positions, api_atoms.positions)
    assert cli_atoms.get_potential_energy() == pytest.approx(api_atoms.get_potential_energy())
    np.testing.assert_allclose(cli_atoms.get_forces(), api_atoms.get_forces())
    assert cli_manifest["stages"] == api_manifest["stages"]
    assert cli_manifest["artifacts"] == api_manifest["artifacts"]
    assert cli_manifest["backend_identity"] == api_manifest["backend_identity"]
    assert cli_result.artifacts == ()
    assert api_result.artifacts == tuple(api_manifest["artifacts"])


def test_api_dispatch_does_not_import_cli_module(tmp_path):
    """Fresh API dispatch resolves the execution layer without loading CLI code."""
    script = """
import sys
from atst_tools.api import RunOptions
from atst_tools.api.services import _dispatch_normalized
from atst_tools.workflows import execution

class FakeWorkflow:
    def __init__(self, config, calculator_name, calculation):
        pass
    def run(self):
        return 'executed'

execution.RelaxWorkflow = FakeWorkflow
assert _dispatch_normalized(
    {'calculation': {'type': 'relax'}, 'calculator': {'name': 'abacus'}},
    RunOptions(),
) == 'executed'
assert 'atst_tools.scripts.main' not in sys.modules
"""
    completed = subprocess.run(
        [sys.executable, "-c", script],
        cwd=tmp_path,
        env=child_test_environment(Path(__file__).resolve().parents[2] / "src"),
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr


def test_cli_and_api_keep_their_error_surfaces(monkeypatch, tmp_path):
    """The same missing input is typed for API and unwrapped for CLI users."""
    from atst_tools.api import RunOptions, run_workflow
    from atst_tools.api.models import WorkflowExecutionError
    from atst_tools.scripts import main as cli_main

    _write_relax_fixture(tmp_path / "missing")
    monkeypatch.chdir(tmp_path / "missing")
    Path("initial.traj").unlink()
    with pytest.raises(WorkflowExecutionError) as api_error:
        run_workflow("config.yaml", RunOptions(world=FakeWorld()))
    with pytest.raises(FileNotFoundError) as cli_error:
        cli_main.main(["config.yaml"])
    assert isinstance(api_error.value.__cause__, FileNotFoundError)
    assert str(api_error.value.__cause__) == str(cli_error.value)


@pytest.mark.parametrize(
    ("configured_vector", "vector_file"),
    [("chosen.npy", "chosen.npy"), (None, "displacement_vector.npy")],
)
def test_shared_dimer_loads_displacement_vector(
    monkeypatch, tmp_path, configured_vector, vector_file
):
    """Explicit and fallback Dimer vectors reach the scientific runner intact."""
    from atst_tools.workflows import execution

    vector = np.array([[0.2, -0.3, 0.4]])
    np.save(tmp_path / vector_file, vector)
    monkeypatch.chdir(tmp_path)
    atoms = Atoms("H", positions=[[0.0, 0.0, 0.0]])
    monkeypatch.setattr(execution, "read_structure", lambda path: atoms)
    seen = {}

    class FakeDimer:
        def __init__(self, **kwargs):
            seen.update(kwargs)

        def run(self, **kwargs):
            seen["run"] = kwargs

    monkeypatch.setattr(execution, "AbacusDimer", FakeDimer)
    result = execution.run_dimer(
        {},
        "abacus",
        {
            "type": "dimer",
            "init_structure": "initial.traj",
            "init_eigenmode_method": "displacement",
            "displacement_vector": configured_vector,
        },
    )

    assert result is atoms
    np.testing.assert_array_equal(seen["displacement_vector"], vector)
    assert "fmax" in seen["run"]
