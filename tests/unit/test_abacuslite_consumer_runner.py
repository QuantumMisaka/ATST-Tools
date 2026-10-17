"""Contract checks for the installed ATST abacuslite parity runner."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


SCRIPT = Path(__file__).resolve().parents[2] / "scripts/run_abacuslite_consumer.py"


def _runner():
    spec = importlib.util.spec_from_file_location("run_abacuslite_consumer", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _report(backend: str) -> dict:
    return {
        "schema": "atst-abacuslite-consumer-v1",
        "backend": backend,
        "packages": {
            "atst": {"version": "2.2.9"}, "python": "3.12.0",
            "atst_files_sha256": "a" * 64, "atst_wheel_sha256": "b" * 64,
            "distributions": {"atst-tools": "2.2.9", "ase": "3.28.0"},
        },
        "tests": {
            "collected": 5,
            "passed": 4,
            "skipped": 1,
            "failed": 0,
            "errors": 0,
            "skip_reasons": {"optional MPI launcher unavailable": 1},
            "testcases": {
                "suite::test_a": {"outcome": "passed"},
                "suite::test_b": {"outcome": "passed"},
                "suite::test_c": {"outcome": "passed"},
                "suite::test_d": {"outcome": "passed"},
                "suite::test_e": {"outcome": "skipped", "reason": "optional MPI launcher unavailable"},
            },
        },
        "facts": {
            "scf": {"energy_eV": -10.0, "forces_eV_per_ang": [[0.5, 0.0, 0.0]]},
            "md": {"frames": 2, "cell_angstrom": [[3.0, 0.0, 0.0]]},
        },
    }


def test_runner_compares_same_collection_and_numeric_facts():
    runner = _runner()
    assert runner.compare_reports(_report("vendored"), _report("external")) == {}


@pytest.mark.parametrize("mutation", [
    lambda report: report["facts"]["scf"].update(energy_eV=-9.99),
    lambda report: report["tests"].update(collected=4),
    lambda report: report["tests"].update(failed=1),
    lambda report: report["tests"]["skip_reasons"].update({"unexpected skip": 1}),
    lambda report: report["packages"].update(atst_files_sha256="different"),
    lambda report: report["packages"].update(distributions={"ase": "3.29.0"}),
    lambda report: report["tests"]["testcases"].update({"suite::test_d": {"outcome": "skipped", "reason": "unexpected"}}),
    lambda report: report["tests"].update(collected=0, passed=0, skipped=0),
])
def test_runner_rejects_parity_mutations(mutation):
    runner = _runner()
    vendored = _report("vendored")
    external = _report("external")
    mutation(external)
    with pytest.raises(ValueError):
        runner.compare_reports(vendored, external)


def test_runner_rejects_source_shadow(tmp_path):
    runner = _runner()
    site = tmp_path / "venv/lib/python3.12/site-packages"
    checkout = tmp_path / "checkout"
    source = checkout / "src/atst_tools/__init__.py"
    with pytest.raises(RuntimeError, match="site-packages"):
        runner.ensure_installed_origin(source, site, checkout)
    installed = site / "atst_tools/__init__.py"
    assert runner.ensure_installed_origin(installed, site, checkout) == installed.resolve()


def test_runner_rejects_report_without_expected_backend():
    runner = _runner()
    report = _report("vendored")
    report["backend"] = "external"
    with pytest.raises(ValueError, match="vendored"):
        runner.compare_reports(report, _report("external"))


def test_external_consumer_guard_rejects_vendored_import_mutation(tmp_path):
    """A consumer regression that imports the snapshot must fail before use."""
    runner = _runner()
    module = tmp_path / "tests/unit/test_backend_routing.py"
    module.parent.mkdir(parents=True)
    module.write_text(
        "def _exercise_external_consumers(scratch):\n"
        "    from atst_tools.external.ASE_interface.abacuslite.io import generalio\n",
        encoding="utf-8",
    )
    with pytest.raises(RuntimeError, match="Vendored backend import") as error:
        runner._guard_runtime_consumers(tmp_path)
    print(str(error.value))
