"""Tests for normalized constant-potential record extraction."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "extract_constant_potential_record.py"
SAI_RESULT = (
    ROOT
    / "docs/reports/data/ATST_CP_RUNTIME_20260921"
    / "cases/cp-gate-single/constant_potential_results.json"
)


def _load_extractor():
    spec = importlib.util.spec_from_file_location(
        "extract_constant_potential_record",
        SCRIPT,
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_extracts_archived_sai_point_without_hand_copying():
    extractor = _load_extractor()
    results = json.loads(SAI_RESULT.read_text(encoding="utf-8"))

    record = extractor.extract_atst_point(results)

    assert record["nelec"] == pytest.approx(217.0)
    assert record["efermi_eV"] == pytest.approx(9.3259944384)
    assert record["mu_calc_eV"] == pytest.approx(15.045970655701026)
    assert record["energy_eV"] == pytest.approx(-39582.384053687325)
    assert len(record["forces_eV_per_A"]) == 12
    assert record["density_integral_e"] == pytest.approx(217.00000000000415)
    assert record["compensation_derivative_eV"] == pytest.approx(5.719976217301026)
    assert record["iterations"] == 1
    assert record["converged"] is True


def test_rejects_missing_or_out_of_range_points():
    extractor = _load_extractor()

    with pytest.raises(ValueError, match="contain no points"):
        extractor.extract_atst_point({"points": []})

    with pytest.raises(ValueError, match="outside"):
        extractor.extract_atst_point({"points": [{"nelec": 1.0}]}, point_index=1)
