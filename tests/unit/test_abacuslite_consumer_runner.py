"""Contract checks for the external-only installed-wheel consumer runner."""

from __future__ import annotations

from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import sys
import sysconfig
from types import ModuleType, SimpleNamespace

import pytest


SCRIPT = Path(__file__).resolve().parents[2] / "scripts/run_abacuslite_consumer.py"
REQUIRED_CASES = (
    "tests.unit.test_backend_routing::test_external_consumers_do_not_import_vendored_backend",
    "tests.unit.test_backend_routing::test_external_frame_selection_fermi_order_and_native_md",
    "tests.unit.test_backend_routing::test_missing_plugin_fails_closed_before_a_vendored_import",
    "tests.unit.test_backend_routing::test_legacy_import_paths_forward_to_the_installed_modules",
    "tests.unit.test_installed_child_provenance::test_child_import_uses_selected_source_and_rejects_installed_shadow",
)


def _runner():
    spec = importlib.util.spec_from_file_location("run_abacuslite_consumer", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _report() -> dict:
    site = "/venv/lib/python3.12/site-packages"
    atst = f"{site}/atst_tools/__init__.py"
    plugin = f"{site}/abacuslite/__init__.py"
    cases = {case: {"outcome": "passed"} for case in REQUIRED_CASES}
    return {
        "schema": "atst-abacuslite-consumer-v2",
        "backend": "external",
        "packages": {
            "atst": {"version": "2.2.9", "path": atst},
            "abacuslite": {"version": "1.1.0", "path": plugin,
                            "files_sha256": "c" * 64, "wheel_sha256": None},
            "python": "3.12.0", "site_packages": site,
            "atst_files_sha256": "a" * 64, "atst_wheel_sha256": "b" * 64,
            "distributions": {"atst-tools": "2.2.9", "ase": "3.28.0"},
        },
        "pytest_cwd": "/checkout",
        "pytest_exit_code": 0,
        "pytest_attestation": {
            phase: {"atst": atst, "abacuslite": plugin, "backend": "external"}
            for phase in ("start", "finish")
        } | {"deselected": []},
        "tests": {
            "collected": len(cases), "passed": len(cases), "skipped": 0,
            "failed": 0, "errors": 0, "deselected": 0,
            "skip_reasons": {}, "testcases": cases,
        },
        "facts": {
            "scf": {"energy_eV": -2898.0, "efermi_eV": -2.04,
                    "forces_eV_per_ang": [[-0.5, -0.01, 0], [0, 0, 0], [0.5, 0.01, 0]]},
            "md": {"frames": 2, "last_energy_eV": -1940.6311063727,
                   "last_forces_eV_per_ang": [
                       [0.1191299488, 0.2678988387, 1.5033416041],
                       [-0.1191299488, -0.2678988387, -1.5033416041]],
                   "cell_angstrom": [[5.34197 * 0.529177, 0, 0], [0, 1, 0], [0, 0, 1]],
                   "pbc": [True, True, True], "symbols": ["As", "Ga"]},
        },
    }


def test_runner_accepts_complete_external_report():
    _runner().verify_report(_report())


def test_runner_clears_inherited_pytest_selection(monkeypatch, tmp_path):
    monkeypatch.setenv("PYTEST_ADDOPTS", "-k test_only_one_case")
    monkeypatch.setenv("PYTEST_PLUGINS", "selective_plugin")
    monkeypatch.setenv("PYTEST_DISABLE_PLUGIN_AUTOLOAD", "1")
    environment = _runner()._pytest_environment(
        tmp_path / "checkout", "external", tmp_path / "attest.json", str(tmp_path / "plugin")
    )
    for key in ("PYTEST_ADDOPTS", "PYTEST_PLUGINS", "PYTEST_DISABLE_PLUGIN_AUTOLOAD"):
        assert key not in environment
    assert environment["ATST_TEST_INSTALLED_WHEEL"] == "1"
    assert environment["ATST_RUNNER_ROOT"] == str(tmp_path / "checkout")


def test_pytest_attestation_starts_fresh_when_report_path_is_reused(monkeypatch, tmp_path):
    runner = _runner()
    site = Path(sysconfig.get_paths()["purelib"])
    atst = ModuleType("atst_tools")
    atst.__file__ = str(site / "atst_tools/__init__.py")
    atst.__path__ = []
    calculators = ModuleType("atst_tools.calculators")
    calculators.__path__ = []
    backend = ModuleType("atst_tools.calculators.abacuslite_backend")
    backend.BACKEND_SOURCE = "external"
    abacuslite = ModuleType("abacuslite")
    abacuslite.__file__ = str(site / "abacuslite/__init__.py")
    for name, module in (("atst_tools", atst),
                         ("atst_tools.calculators", calculators),
                         ("atst_tools.calculators.abacuslite_backend", backend),
                         ("abacuslite", abacuslite)):
        monkeypatch.setitem(sys.modules, name, module)
    attestation = tmp_path / "consumer.pytest.attest.json"
    monkeypatch.setenv("ATST_RUNNER_ROOT", str(tmp_path / "checkout"))
    monkeypatch.setenv("ATST_RUNNER_BACKEND", "external")
    monkeypatch.setenv("ATST_RUNNER_ATTEST", str(attestation))
    hooks = {}
    exec(runner._PYTEST_PROVENANCE_PLUGIN, hooks)

    hooks["pytest_sessionstart"](None)
    hooks["pytest_deselected"]([SimpleNamespace(nodeid="tests/test_sample.py::test_omitted")])
    hooks["pytest_sessionfinish"](None, 0)
    assert json.loads(attestation.read_text())["deselected"] == [
        "tests/test_sample.py::test_omitted"
    ]

    hooks["pytest_sessionstart"](None)
    hooks["pytest_sessionfinish"](None, 0)
    second = json.loads(attestation.read_text())
    assert set(second) == {"start", "finish", "deselected"}
    assert second["deselected"] == []


@pytest.mark.parametrize("mutation", [
    lambda report: report.update(backend="vendored"),
    lambda report: report["facts"]["scf"].update(energy_eV=-2897.0),
    lambda report: report["tests"].update(collected=0, passed=0),
    lambda report: report["tests"].update(failed=1),
    lambda report: report["tests"]["testcases"].pop(REQUIRED_CASES[0]),
    lambda report: report["tests"]["testcases"][REQUIRED_CASES[0]].update(outcome="skipped"),
    lambda report: report["packages"]["abacuslite"].update(path="/checkout/src/abacuslite/__init__.py"),
    lambda report: report["packages"]["abacuslite"].update(files_sha256="missing"),
    lambda report: report["packages"].update(site_packages=""),
    lambda report: report.update(pytest_cwd=""),
    lambda report: report["pytest_attestation"]["finish"].update(backend="vendored"),
    lambda report: report.update(pytest_exit_code=1),
    lambda report: report["tests"].update(deselected=1),
    lambda report: report["pytest_attestation"].update(deselected=["tests/unit/test_omitted.py::test_omitted"]),
    lambda report: report["tests"].pop("deselected"),
])
def test_runner_rejects_external_report_mutations(mutation):
    report = deepcopy(_report())
    mutation(report)
    with pytest.raises((ValueError, RuntimeError)):
        _runner().verify_report(report)


def test_runner_rejects_source_shadow(tmp_path):
    runner = _runner()
    site = tmp_path / "venv/lib/python3.12/site-packages"
    checkout = tmp_path / "checkout"
    source = checkout / "src/atst_tools/__init__.py"
    with pytest.raises(RuntimeError, match="site-packages"):
        runner.ensure_installed_origin(source, site, checkout)
    installed = site / "atst_tools/__init__.py"
    assert runner.ensure_installed_origin(installed, site, checkout) == installed.resolve()


def test_external_consumer_guard_rejects_legacy_runtime_import_mutation(tmp_path):
    """A production regression that imports the old path must fail before use."""
    runner = _runner()
    module = tmp_path / "tests/unit/test_backend_routing.py"
    module.parent.mkdir(parents=True)
    module.write_text(
        "def _exercise_external_consumers(scratch):\n"
        "    from atst_tools.external.ASE_interface.abacuslite.io import generalio\n",
        encoding="utf-8",
    )
    with pytest.raises(RuntimeError, match="Vendored backend import"):
        runner._guard_runtime_consumers(tmp_path)
