"""The ATST CI must exercise the required installed abacuslite backend."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
WORKFLOW = ROOT / ".github/workflows/tests.yml"


def test_installed_consumer_ci_has_one_external_lane_and_verifies_its_report():
    workflow = WORKFLOW.read_text(encoding="utf-8")
    assert "--backend external" in workflow
    assert "--verify-report" in workflow
    assert "--backend vendored" not in workflow
    assert "venv-vendored" not in workflow
    assert "--compare" not in workflow


def test_ci_keeps_minimum_python_ase_matrix_and_no_live_snapshot_drift_job():
    workflow = WORKFLOW.read_text(encoding="utf-8")
    assert 'python-version: ["3.10", "3.12"]' in workflow
    assert 'ase-version: ["3.28.0", "latest"]' in workflow
    assert not (ROOT / ".github/workflows/abacuslite-ase-interface.yml").exists()
