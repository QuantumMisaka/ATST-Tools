"""Tests for the constant-potential result comparator."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "compare_constant_potential_results.py"


def _load_comparator():
    spec = importlib.util.spec_from_file_location("compare_constant_potential_results", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _record() -> dict:
    return {
        "nelec": 216.5,
        "efermi_eV": -3.2,
        "vacuum_level_eV": 1.0,
        "fermishift_eV": -1.0,
        "mu_calc_eV": -4.2,
        "mu_target_eV": -4.2,
        "residual_mu_eV": 0.0,
        "energy_eV": -1000.0,
        "free_energy_eV": -1000.0,
        "forces_eV_per_A": [[0.0, 0.0, 0.0], [0.1, -0.2, 0.3]],
        "density_integral_e": 216.5,
        "compensation_derivative_eV": 0.5,
        "capacitance_e_per_V": 0.02,
        "iterations": 5,
        "converged": True,
    }


def test_identical_records_pass():
    comparator = _load_comparator()
    report = comparator.compare_records(_record(), _record())

    assert report["ok"] is True
    assert report["mismatches"] == []
    assert "forces_eV_per_A" in report["compared"]


def test_scalar_difference_beyond_tolerance_is_reported():
    comparator = _load_comparator()
    candidate = _record()
    candidate["energy_eV"] += 1e-4

    report = comparator.compare_records(_record(), candidate)

    assert report["ok"] is False
    assert any(item["field"] == "energy_eV" for item in report["mismatches"])


def test_force_difference_and_shape_mismatch_are_reported():
    comparator = _load_comparator()
    candidate = _record()
    candidate["forces_eV_per_A"][0][1] = 1e-3

    report = comparator.compare_records(_record(), candidate)
    assert report["ok"] is False
    assert any("forces_eV_per_A[0][1]" == item["field"] for item in report["mismatches"])

    candidate["forces_eV_per_A"] = candidate["forces_eV_per_A"][:1]
    report = comparator.compare_records(_record(), candidate)
    assert any(item["field"] == "forces_eV_per_A.shape" for item in report["mismatches"])


def test_required_field_and_tolerance_override():
    comparator = _load_comparator()
    candidate = _record()
    del candidate["efermi_eV"]

    with pytest.raises(ValueError, match="candidate:efermi_eV"):
        comparator.compare_records(_record(), candidate)

    candidate = _record()
    candidate["energy_eV"] += 1e-4
    report = comparator.compare_records(
        _record(), candidate, {"energy_eV": 1e-3}
    )
    assert report["ok"] is True


def test_iterations_and_convergence_require_exact_agreement():
    comparator = _load_comparator()
    candidate = _record()
    candidate["iterations"] = 6
    candidate["converged"] = False

    report = comparator.compare_records(_record(), candidate)

    assert report["ok"] is False
    fields = {item["field"] for item in report["mismatches"]}
    assert "converged" in fields
    informational = {item["field"] for item in report["informational"]}
    assert "iterations" in informational


def test_capacitance_is_informational():
    comparator = _load_comparator()
    candidate = _record()
    candidate["capacitance_e_per_V"] = 0.03

    report = comparator.compare_records(_record(), candidate)

    assert report["ok"] is True
    assert any(item["field"] == "capacitance_e_per_V" for item in report["informational"])
