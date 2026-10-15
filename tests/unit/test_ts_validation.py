"""Scope and evidence tests for transition-state screening."""

import json

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.constraints import FixAtoms
from ase.io import write
from ase.vibrations import Vibrations

from atst_tools.api.validation import validate_ts
from atst_tools.utils.ts_validation import build_ts_validation_summary, frozen_structure_fact


GEOMETRY = {
    "numbers": [1, 1, 1],
    "positions": [[0, 0, 0], [1, 0, 0], [2, 0, 0]],
    "cell": [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
    "pbc": [False, False, False],
}


def test_legacy_screen_pass_does_not_certify_ts_without_stationarity():
    summary = build_ts_validation_summary(
        {"imaginary_frequencies": [120.0, 0.0], "indices": [0]},
    )

    assert summary["status"] == "pass"
    assert summary["checks"]["n_imaginary_modes"]["pass"] is True
    assert summary["checks"]["fmax"]["pass"] is None
    assert summary["stationarity"]["status"] == "unknown"
    assert summary["requested_checks_status"] == "unknown"
    assert summary["scope"]["active_indices"] == [0]
    assert summary["scope"]["covers_all_unconstrained_atoms"] is None


def test_frozen_forces_are_checked_only_on_declared_active_scope():
    summary = build_ts_validation_summary(
        {"imaginary_frequencies": [150.0], "indices": [1, 2], "n_atoms": 3, "geometry": GEOMETRY},
        structure={"n_atoms": 3, "geometry": GEOMETRY, "fixed_indices": [0], "forces": [[9.0, 0, 0], [0.02, 0, 0], [0, 0.03, 0]]},
    )

    assert summary["scope"]["kind"] == "constrained"
    assert summary["scope"]["active_indices"] == [1, 2]
    assert summary["scope"]["covers_all_unconstrained_atoms"] is True
    assert summary["stationarity"]["status"] == "pass"
    assert summary["stationarity"]["fmax"] == pytest.approx(0.03)
    assert summary["requested_checks_status"] == "pass"
    assert summary["scope"]["full_system"] is False


def test_partial_vibration_and_missing_constraints_keep_comprehensive_check_unknown():
    summary = build_ts_validation_summary(
        {"imaginary_frequencies": [150.0], "indices": [1], "n_atoms": 3, "geometry": GEOMETRY},
        structure={"n_atoms": 3, "geometry": GEOMETRY, "forces": [[0, 0, 0], [0.02, 0, 0], [0, 0, 0]]},
    )

    assert summary["scope"]["kind"] == "subset"
    assert summary["scope"]["covers_all_unconstrained_atoms"] is None
    assert summary["stationarity"]["status"] == "pass"
    assert summary["requested_checks_status"] == "unknown"


def test_missing_vibration_indices_never_imply_full_system_scope():
    summary = build_ts_validation_summary(
        {"imaginary_frequencies": [120.0], "n_atoms": 3, "geometry": GEOMETRY},
        structure={
            "n_atoms": 3,
            "geometry": GEOMETRY,
            "fixed_indices": [],
            "forces": [[0.01, 0, 0]] * 3,
        },
    )

    assert summary["scope"]["active_indices"] is None
    assert summary["scope"]["full_system"] is False
    assert summary["scope"]["covers_all_unconstrained_atoms"] is None
    assert summary["stationarity"]["status"] == "unknown"
    assert summary["requested_checks_status"] == "unknown"


def test_explicit_full_selection_is_order_independent():
    geometry = {
        "numbers": [1, 1],
        "positions": [[0, 0, 0], [1, 0, 0]],
        "cell": [[0, 0, 0]] * 3,
        "pbc": [False] * 3,
    }
    summary = build_ts_validation_summary(
        {"imaginary_frequencies": [120.0], "indices": [1, 0], "n_atoms": 2, "geometry": geometry},
        structure={"n_atoms": 2, "geometry": geometry, "fixed_indices": [], "forces": [[0.01, 0, 0]] * 2},
    )

    assert summary["scope"]["active_indices"] == [1, 0]
    assert summary["scope"]["full_system"] is True
    assert summary["requested_checks_status"] == "pass"


def test_explicit_selection_overlapping_fixed_atoms_stays_screenable_but_uncertain():
    atoms = _constrained_frozen_atoms()
    fact = frozen_structure_fact(atoms)
    summary = build_ts_validation_summary(
        {
            "imaginary_frequencies": [120.0],
            "indices": [0, 1],
            "n_atoms": 2,
            "geometry": fact["geometry"],
            "fixed_indices": [0],
        },
        structure=fact,
    )

    assert summary["status"] == "pass"
    assert summary["scope"]["active_indices"] == [0, 1]
    assert summary["scope"]["fixed_indices_in_selection"] == [0]
    assert summary["scope"]["full_system"] is False
    assert summary["stationarity"]["status"] == "unknown"
    assert summary["requested_checks_status"] == "unknown"


@pytest.mark.parametrize(
    "results, message",
    [
        ({"imaginary_frequencies": [float("nan")]}, "finite"),
        ({"imaginary_frequencies": "120"}, "sequence"),
        ({"imaginary_frequencies": [100], "indices": [3], "n_atoms": 2}, "out of range"),
        ({"imaginary_frequencies": [100], "indices": [True]}, "integer"),
    ],
)
def test_malformed_vibration_facts_are_diagnostic(results, message):
    with pytest.raises((TypeError, ValueError), match=message):
        build_ts_validation_summary(results)


def test_validate_ts_reads_existing_json_and_never_uses_live_calculator(tmp_path):
    vibration_path = tmp_path / "vibration.json"
    structure_path = tmp_path / "structure.json"
    output_path = tmp_path / "validation.json"
    geometry = {
        "numbers": [1], "positions": [[0, 0, 0]],
        "cell": [[0, 0, 0], [0, 0, 0], [0, 0, 0]], "pbc": [False, False, False],
    }
    vibration_path.write_text(json.dumps({"imaginary_frequencies": [120.0], "indices": [0], "n_atoms": 1, "geometry": geometry}))
    structure_path.write_text(json.dumps({"n_atoms": 1, "geometry": geometry, "fixed_indices": [], "forces": [[0.01, 0, 0]]}))

    summary = validate_ts(vibration_path, structure=structure_path, output=output_path)

    assert summary["requested_checks_status"] == "pass"
    assert summary["scope"]["full_system"] is True
    assert json.loads(output_path.read_text()) == summary


def test_missing_imaginary_modes_never_passes_requested_checks():
    summary = build_ts_validation_summary({"indices": [0], "n_atoms": 1})
    assert summary["status"] == "fail"
    assert summary["checks"]["n_imaginary_modes"]["pass"] is None
    assert summary["requested_checks_status"] == "unknown"


def test_unpaired_force_fact_does_not_claim_stationarity():
    summary = build_ts_validation_summary(
        {"imaginary_frequencies": [120.0], "n_atoms": 1},
        structure={"n_atoms": 1, "forces": [[0, 0, 0]], "fixed_indices": []},
    )
    assert summary["stationarity"]["status"] == "unknown"
    assert summary["requested_checks_status"] == "unknown"


def test_ase_trajectory_uses_only_paired_single_point_forces(tmp_path):
    atoms = Atoms("H", positions=[[0, 0, 0]])
    atoms.calc = SinglePointCalculator(atoms, energy=0.0, forces=[[0.01, 0, 0]])
    trajectory = tmp_path / "frozen.traj"
    write(trajectory, atoms)
    vibration = tmp_path / "vibration.json"
    fact = frozen_structure_fact(atoms)
    vibration.write_text(json.dumps({
        "imaginary_frequencies": [120.0], "n_atoms": 1, "indices": [0],
        "fixed_indices": [], "geometry": fact["geometry"],
    }))

    summary = validate_ts(vibration, structure=trajectory)

    assert summary["stationarity"]["status"] == "pass"
    assert summary["requested_checks_status"] == "pass"


def test_mismatched_geometry_rejects_force_evidence():
    shifted = dict(GEOMETRY, positions=[[1, 0, 0], [1, 0, 0], [2, 0, 0]])
    with pytest.raises(ValueError, match="geometry does not match"):
        build_ts_validation_summary(
            {"imaginary_frequencies": [120.0], "n_atoms": 3, "geometry": GEOMETRY},
            structure={"n_atoms": 3, "geometry": shifted, "forces": [[0, 0, 0]] * 3},
        )


def _constrained_frozen_atoms():
    atoms = Atoms("H2", positions=[[0, 0, 0], [1, 0, 0]])
    atoms.set_constraint(FixAtoms(indices=[0]))
    atoms.calc = SinglePointCalculator(atoms, energy=0.0, forces=[[9.0, 0, 0], [0.01, 0, 0]])
    return atoms


class _NoRunVibrations(Vibrations):
    """Keep ASE's scope resolution while avoiding displacement calculations."""

    def run(self):
        return None

    def summary(self):
        return None

    def get_energies(self):
        return np.array([0.1])

    def get_frequencies(self):
        return np.array([120j])

    def get_zero_point_energy(self):
        return 0.05


def test_vibration_workflow_persists_ase_resolved_constraint_scope(monkeypatch, tmp_path):
    from atst_tools.workflows import vibration

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(vibration, "read_structure", lambda path: _constrained_frozen_atoms())
    monkeypatch.setattr(vibration.CalculatorFactory, "get_calculator", lambda *args, **kwargs: None)
    monkeypatch.setattr(vibration, "Vibrations", _NoRunVibrations)
    monkeypatch.setattr(vibration, "compute_vibration_thermochemistry", lambda *args: {"model": "harmonic", "temperature": 300.0})
    (tmp_path / "ts.traj").touch()

    workflow = vibration.VibrationWorkflow(
        {"calculator": {"name": "dp", "dp": {"model": "unused"}}},
        "dp",
        {"type": "vibration", "init_structure": "ts.traj"},
    )
    workflow.run()

    results = json.loads((tmp_path / "vibration_results.json").read_text(encoding="utf-8"))
    validation = json.loads((tmp_path / "ts_validation.json").read_text(encoding="utf-8"))
    assert results["indices"] == [1]
    assert validation["scope"]["active_indices"] == [1]
    assert validation["stationarity"]["fmax"] == pytest.approx(0.01)
    assert validation["requested_checks_status"] == "pass"


def test_vibration_workflow_accepts_explicit_fixed_atom_selection(monkeypatch, tmp_path):
    from atst_tools.workflows import vibration

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(vibration, "read_structure", lambda path: _constrained_frozen_atoms())
    monkeypatch.setattr(vibration.CalculatorFactory, "get_calculator", lambda *args, **kwargs: None)
    monkeypatch.setattr(vibration, "Vibrations", _NoRunVibrations)
    monkeypatch.setattr(vibration, "compute_vibration_thermochemistry", lambda *args: {"model": "harmonic", "temperature": 300.0})
    (tmp_path / "ts.traj").touch()
    workflow = vibration.VibrationWorkflow(
        {"calculator": {"name": "dp", "dp": {"model": "unused"}}},
        "dp",
        {"type": "vibration", "init_structure": "ts.traj", "indices": [0, 1]},
    )

    workflow.run()

    results = json.loads((tmp_path / "vibration_results.json").read_text(encoding="utf-8"))
    validation = json.loads((tmp_path / "ts_validation.json").read_text(encoding="utf-8"))
    assert results["indices"] == [0, 1]
    assert validation["scope"]["fixed_indices_in_selection"] == [0]
    assert validation["stationarity"]["status"] == "unknown"
    assert validation["requested_checks_status"] == "unknown"


def test_d2s_workflow_persists_ase_resolved_constraint_scope(monkeypatch, tmp_path):
    from atst_tools.workflows import d2s

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(d2s, "Vibrations", _NoRunVibrations)
    monkeypatch.setattr(d2s, "compute_vibration_thermochemistry", lambda *args: {"model": "harmonic", "temperature": 300.0})
    monkeypatch.setattr(d2s.D2SWorkflow, "_get_calc", lambda *args: None)
    workflow = d2s.D2SWorkflow(
        {"calculator": {"name": "dp", "dp": {"model": "unused"}}},
        "dp",
        {"type": "d2s", "method": "dimer", "vibration": {"enabled": True, "indices": "all"}},
    )
    frozen_path = tmp_path / "frozen.traj"
    write(frozen_path, _constrained_frozen_atoms())
    workflow.run_vibration([], _constrained_frozen_atoms(), str(frozen_path))

    results = json.loads((tmp_path / "d2s_vibration_results.json").read_text(encoding="utf-8"))
    validation = json.loads((tmp_path / "d2s_ts_validation.json").read_text(encoding="utf-8"))
    assert results["indices"] == [1]
    assert validation["scope"]["active_indices"] == [1]
    assert validation["stationarity"]["fmax"] == pytest.approx(0.01)
    assert validation["requested_checks_status"] == "pass"
