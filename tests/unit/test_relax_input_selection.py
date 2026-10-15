"""Regression tests for Relax's initial-structure selection."""

from __future__ import annotations

import pytest
import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator
from ase.io import read, write

from atst_tools.workflows import relax


class PositionCalculator(Calculator):
    """Record the geometry evaluated by a controlled ASE calculator."""

    implemented_properties = ("energy", "forces")

    def __init__(self):
        super().__init__()
        self.positions = []

    def calculate(self, atoms=None, properties=("energy",), system_changes=()):
        super().calculate(atoms, properties, system_changes)
        self.positions.append(float(atoms.positions[0, 0]))
        self.results = {"energy": 0.0, "forces": np.zeros((len(atoms), 3))}


def _run_relax(monkeypatch, calculator, *, restart=False):
    """Run Relax with its real ASE optimizer and a controlled calculator."""
    monkeypatch.setattr(
        relax.CalculatorFactory,
        "get_calculator",
        lambda name, config, **kwargs: calculator,
    )
    return relax.RelaxWorkflow(
        {},
        "controlled",
        {
            "type": "relax",
            "init_structure": "configured.traj",
            "trajectory": "relax.traj",
            "logfile": "relax.log",
            "optimizer": "BFGS",
            "fmax": 0.1,
            "max_steps": 1,
            "restart": restart,
            "artifact_manifest": None,
        },
    ).run()


@pytest.mark.parametrize(
    ("available", "expected"),
    [
        ({"configured.traj": 0.25, "init.traj": 0.75}, 0.25),
        ({"init.traj": 0.75}, 0.75),
    ],
)
def test_relax_uses_configured_structure_or_cwd_fallback(
    monkeypatch, tmp_path, available, expected
):
    """A missing configured file falls back to the current directory's init.traj."""
    monkeypatch.chdir(tmp_path)
    for filename, position in available.items():
        write(filename, Atoms("H", positions=[[position, 0.0, 0.0]]))
    calculator = PositionCalculator()

    _run_relax(monkeypatch, calculator)

    assert calculator.positions == [expected]
    assert read("final_relaxed.traj").positions[0, 0] == expected


def test_relax_reports_missing_initial_structure(monkeypatch, tmp_path):
    """No configured or fallback input fails before creating a calculator."""
    monkeypatch.chdir(tmp_path)
    calculator = PositionCalculator()

    with pytest.raises(FileNotFoundError, match="configured.traj"):
        _run_relax(monkeypatch, calculator)

    assert calculator.positions == []


def test_relax_restart_uses_last_trajectory_frame(monkeypatch, tmp_path):
    """A restart uses the last saved frame despite an available input file."""
    monkeypatch.chdir(tmp_path)
    write("configured.traj", Atoms("H", positions=[[0.25, 0.0, 0.0]]))
    write(
        "relax.traj",
        [
            Atoms("H", positions=[[0.5, 0.0, 0.0]]),
            Atoms("H", positions=[[0.75, 0.0, 0.0]]),
        ],
    )
    calculator = PositionCalculator()

    _run_relax(monkeypatch, calculator, restart=True)

    assert calculator.positions == [0.75]
    assert read("final_relaxed.traj").positions[0, 0] == 0.75
