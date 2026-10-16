"""Regression tests for IRC trajectory paths and cached descent results."""

import json

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read, write
from ase.optimize import FIRE

from atst_tools.workflows import irc


class CountingHarmonic(Calculator):
    """Provide a deterministic minimum while counting backend evaluations."""

    implemented_properties = ["energy", "forces"]

    def __init__(self):
        super().__init__()
        self.calls = 0

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self.calls += 1
        positions = atoms.get_positions()
        self.results["energy"] = float(np.sum(positions**2))
        self.results["forces"] = -2 * positions


@pytest.mark.parametrize("normalized_setting", ["omitted", None, "custom-normalized.traj"])
def test_sella_both_writes_normalized_trajectory_for_default_null_and_custom(
    tmp_path, monkeypatch, normalized_setting
):
    """Both completed directions yield a usable normalized artifact path."""
    import sella

    class CompletedIRC:
        def __init__(self, atoms, trajectory, **kwargs):
            self.atoms = atoms
            self.trajectory = trajectory
            self.nsteps = 0

        def run(self, fmax, steps, direction):
            self.atoms.set_positions([[0.2 if direction == "forward" else -0.2, 0, 0]])
            self.atoms.get_potential_energy()
            self.atoms.get_forces()
            self.trajectory.write(self.atoms)
            self.nsteps += 1
            return True

    monkeypatch.chdir(tmp_path)
    write("start.traj", Atoms("H", positions=[[0, 0, 0]]))
    monkeypatch.setattr(sella, "IRC", CompletedIRC)
    monkeypatch.setattr(
        irc.CalculatorFactory, "get_calculator", lambda *args, **kwargs: CountingHarmonic()
    )
    calculation = {"type": "irc", "init_structure": "start.traj", "trajectory": "irc.traj"}
    if normalized_setting != "omitted":
        calculation["normalized_trajectory"] = normalized_setting

    workflow = irc.IRCWorkflow({}, "dp", calculation)
    workflow.run()

    expected = (
        normalized_setting if normalized_setting == "custom-normalized.traj" else "norm_irc.traj"
    )
    assert workflow.normalized_traj_file == expected
    normalized = read(expected, index=":")
    assert len(normalized) == 2
    assert all(np.isfinite(frame.get_potential_energy()) for frame in normalized)
    manifest = json.loads((tmp_path / "atst_artifacts.json").read_text())
    assert {artifact["path"] for artifact in manifest["artifacts"]} == {"irc.traj", expected}


def test_descent_preserves_cached_results_without_save_time_evaluation(tmp_path, monkeypatch):
    """Every written endpoint has its own frozen energy and forces."""
    monkeypatch.chdir(tmp_path)
    write("start.traj", Atoms("H", positions=[[0, 0, 0]]))
    np.save("mode.npy", np.array([[1.0, 0.0, 0.0]]))
    calculators = []
    completion_counts = []

    class ObservedFIRE(FIRE):
        def run(self, *args, **kwargs):
            result = super().run(*args, **kwargs)
            completion_counts.append(sum(calculator.calls for calculator in calculators))
            return result

    def make_calculator(*args, **kwargs):
        calculator = CountingHarmonic()
        calculators.append(calculator)
        return calculator

    monkeypatch.setattr(irc.CalculatorFactory, "get_calculator", make_calculator)
    monkeypatch.setattr(irc, "FIRE", ObservedFIRE)
    original_write = irc.write
    save_counts = []

    def observing_write(*args, **kwargs):
        before = sum(calculator.calls for calculator in calculators)
        result = original_write(*args, **kwargs)
        save_counts.append((before, sum(calculator.calls for calculator in calculators)))
        return result

    monkeypatch.setattr(irc, "write", observing_write)
    workflow = irc.IRCWorkflow(
        {}, "dp", {"type": "irc", "backend": "descent", "init_structure": "start.traj",
                    "mode_vector": "mode.npy", "trajectory": "irc.traj",
                    "normalized_trajectory": "custom-normalized.traj", "direction": "both",
                    "descent_delta": 0.2, "max_steps": 1, "fmax": 0.01},
    )
    returned = workflow.run()

    assert len(calculators) == 2
    assert sum(calculator.calls for calculator in calculators) == completion_counts[-1]
    assert len(save_counts) == 2 and all(before == after for before, after in save_counts)
    main = read("irc.traj", index=":")
    normalized = read("custom-normalized.traj", index=":")
    assert len(returned) == len(main) == len(normalized) == 2
    for index, direction in enumerate(("forward", "reverse")):
        branch = read(f"irc_{direction}.traj", index=-1)
        for frame in (returned[index], main[index], normalized[index], branch):
            assert isinstance(frame.calc, SinglePointCalculator)
            np.testing.assert_allclose(frame.positions, branch.positions)
            np.testing.assert_allclose(frame.get_potential_energy(), branch.get_potential_energy())
            np.testing.assert_allclose(frame.get_forces(), branch.get_forces())
        assert returned[index].calc is not branch.calc


@pytest.mark.parametrize("recorded_position", [0.2, 0.3])
def test_descent_snapshot_reuses_only_matching_recorded_endpoint(tmp_path, recorded_position):
    """Recover a missing cached property only from the same stored geometry."""
    branch = Atoms("H", positions=[[0.2, 0, 0]])
    branch.calc = SinglePointCalculator(branch, forces=np.array([[-0.4, 0, 0]]))
    recorded = Atoms("H", positions=[[recorded_position, 0, 0]])
    recorded.calc = SinglePointCalculator(
        recorded, energy=recorded_position**2, forces=np.array([[-0.4, 0, 0]])
    )
    trajectory = tmp_path / "branch.traj"
    write(trajectory, recorded)

    snapshot = irc.IRCWorkflow._freeze_descent_frame(branch, str(trajectory))

    assert snapshot.calc is not branch.calc
    if recorded_position == 0.2:
        np.testing.assert_allclose(snapshot.get_potential_energy(), 0.04)
        np.testing.assert_allclose(snapshot.get_forces(), [[-0.4, 0, 0]])
    else:
        assert "energy" not in snapshot.calc.results
