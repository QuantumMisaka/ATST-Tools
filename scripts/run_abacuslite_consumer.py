#!/usr/bin/env python3
"""Run ATST's installed-wheel suite against its required external abacuslite.

The caller installs both wheels and invokes this script from outside their
checkouts. The suite uses curated raw fixtures; no engine or network is used.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import importlib
from importlib import metadata
from importlib.util import find_spec
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import sysconfig
import tempfile
import xml.etree.ElementTree as ET


SCHEMA = "atst-abacuslite-consumer-v2"
_REQUIRED_CASES = (
    "tests.unit.test_backend_routing::test_external_consumers_do_not_import_vendored_backend",
    "tests.unit.test_backend_routing::test_external_frame_selection_fermi_order_and_native_md",
    "tests.unit.test_backend_routing::test_missing_plugin_fails_closed_before_a_vendored_import",
    "tests.unit.test_backend_routing::test_legacy_import_paths_forward_to_the_installed_modules",
    "tests.unit.test_installed_child_provenance::test_child_import_uses_selected_source_and_rejects_installed_shadow",
)
# Literal lattice from multiframe_md_latest/MD_dump: LATTICE_VECTORS times
# LATTICE_CONSTANT (Angstrom). Keep this oracle independent of plugin readers.
_MD_CELL_ANGSTROM = [
    [5.34197 * 0.529177, 5.34197 * 0.529177, 0.0],
    [0.0, 5.34197 * 0.529177, 5.34197 * 0.529177],
    [5.34197 * 0.529177, 0.0, 5.34197 * 0.529177],
]


def ensure_installed_origin(path: str | Path, site_packages: Path, checkout: Path) -> Path:
    """Require a wheel import from this interpreter's site-packages directory."""
    resolved = Path(path).resolve()
    site_packages = site_packages.resolve()
    checkout = checkout.resolve()
    if resolved.is_relative_to(checkout) or not resolved.is_relative_to(site_packages):
        raise RuntimeError(f"Expected installed site-packages import, got {resolved}")
    return resolved


def _compare_values(left, right, name: str) -> None:
    if isinstance(left, dict) and isinstance(right, dict):
        if left.keys() != right.keys():
            raise ValueError(f"Different fact keys at {name}: {left.keys()} != {right.keys()}")
        for key in left:
            _compare_values(left[key], right[key], f"{name}.{key}")
    elif isinstance(left, list) and isinstance(right, list):
        if len(left) != len(right):
            raise ValueError(f"Different fact lengths at {name}: {len(left)} != {len(right)}")
        for index, (a, b) in enumerate(zip(left, right)):
            _compare_values(a, b, f"{name}[{index}]")
    elif isinstance(left, bool) or isinstance(right, bool):
        if left is not right:
            raise ValueError(f"Different Boolean fact at {name}: {left!r} != {right!r}")
    elif isinstance(left, (int, float)) and isinstance(right, (int, float)):
        if not (math.isfinite(left) and math.isfinite(right)
                and math.isclose(left, right, rel_tol=0, abs_tol=1e-8)):
            raise ValueError(f"Different numeric fact at {name}: {left!r} != {right!r}")
    elif left != right:
        raise ValueError(f"Different fact at {name}: {left!r} != {right!r}")


def verify_report(report: dict) -> None:
    """Fail closed on incomplete external-suite evidence and source provenance."""
    if report.get("schema") != SCHEMA or report.get("backend") != "external":
        raise ValueError(f"Expected external {SCHEMA} report")
    packages = report.get("packages", {})
    plugin = packages.get("abacuslite") or {}
    if not report.get("pytest_cwd") or not packages.get("site_packages"):
        raise ValueError("Missing pytest checkout or site-packages provenance")
    root = Path(report["pytest_cwd"])
    site = Path(packages["site_packages"])
    if not root.is_absolute() or not site.is_absolute() or site == Path("/"):
        raise ValueError("Invalid pytest checkout or site-packages provenance")
    for name, item in (("ATST", packages.get("atst") or {}),
                       ("abacuslite", plugin)):
        if not item.get("version") or not item.get("path"):
            raise ValueError(f"Missing {name} installed distribution provenance")
        ensure_installed_origin(item["path"], site, root)
    if not packages.get("python") or not packages.get("distributions"):
        raise ValueError("Missing Python or installed dependency provenance")
    for name in ("atst_files_sha256", "atst_wheel_sha256"):
        if len(packages.get(name) or "") != 64:
            raise ValueError(f"Missing {name} provenance")
    if len(plugin.get("files_sha256") or "") != 64:
        raise ValueError("Missing abacuslite installed-file provenance")
    plugin_wheel = plugin.get("wheel_sha256")
    if plugin_wheel is not None and len(plugin_wheel) != 64:
        raise ValueError("Malformed abacuslite wheel provenance")
    attest = report.get("pytest_attestation") or {}
    if attest.get("deselected") != []:
        raise ValueError("External offline suite deselected tests")
    for phase in ("start", "finish"):
        phase_record = attest.get(phase) or {}
        if (phase_record.get("backend") != "external"
                or phase_record.get("atst") != packages["atst"]["path"]
                or phase_record.get("abacuslite") != plugin["path"]):
            raise ValueError(f"Missing or inconsistent pytest {phase} wheel attestation")
    tests = report.get("tests") or {}
    if tests.get("deselected") != 0:
        raise ValueError("External offline suite deselected tests")
    collected = tests.get("collected", 0)
    cases = tests.get("testcases") or {}
    if (collected <= 0 or len(cases) != collected
            or tests.get("passed", 0) + tests.get("skipped", 0) != collected
            or tests.get("failed") or tests.get("errors")
            or report.get("pytest_exit_code") != 0):
        raise ValueError("External offline suite is missing, incomplete or failed")
    if sum(tests.get("skip_reasons", {}).values()) != tests.get("skipped"):
        raise ValueError("External offline skip reasons are incomplete")
    if any("no-plugin environment" in reason or "external-package lane" in reason
           for reason in tests.get("skip_reasons", {})):
        raise ValueError("External backend tests were unexpectedly skipped")
    for case in _REQUIRED_CASES:
        if cases.get(case) != {"outcome": "passed"}:
            raise ValueError(f"Required external consumer case did not pass: {case}")
    facts = report.get("facts") or {}
    expected = {
        "scf": {"energy_eV": -2898.0, "efermi_eV": -2.04,
                "forces_eV_per_ang": [[-0.5, -0.01, 0], [0, 0, 0], [0.5, 0.01, 0]]},
        "md": {"frames": 2, "last_energy_eV": -1940.6311063727,
               "last_forces_eV_per_ang": [[0.1191299488, 0.2678988387, 1.5033416041],
                                             [-0.1191299488, -0.2678988387, -1.5033416041]],
               "cell_angstrom": _MD_CELL_ANGSTROM,
               "pbc": [True, True, True], "symbols": ["As", "Ga"]},
    }
    for section, values in expected.items():
        for key, value in values.items():
            _compare_values(facts.get(section, {}).get(key), value, f"facts.{section}.{key}")


def _installed_file_digest(distribution: metadata.Distribution, prefix: str) -> str:
    """Fingerprint installed production files, independently of wheel location."""
    files = sorted(str(file) for file in (distribution.files or ())
                   if str(file).startswith(prefix) and str(file).endswith(".py"))
    if not files:
        raise RuntimeError(f"No installed production files found under {prefix}")
    digest = hashlib.sha256()
    for name in files:
        digest.update(name.encode() + b"\0")
        digest.update(Path(distribution.locate_file(name)).read_bytes())
        digest.update(b"\0")
    return digest.hexdigest()


def _installed_wheel_digest(distribution: metadata.Distribution) -> str | None:
    """Record the exact wheel hash when pip preserved its direct URL metadata."""
    source = distribution.read_text("direct_url.json")
    if not source:
        return None
    archive_info = json.loads(source).get("archive_info", {})
    hashes = archive_info.get("hashes", {})
    if "sha256" in hashes:
        return hashes["sha256"]
    legacy = archive_info.get("hash", "")
    return legacy.removeprefix("sha256=") if legacy.startswith("sha256=") else None


def _provenance(root: Path, backend: str) -> dict:
    if backend != "external":
        raise RuntimeError("ATST requires the external abacuslite backend")
    if Path.cwd().resolve().is_relative_to(root):
        raise RuntimeError("Run installed consumer checks from outside the ATST checkout")
    site_packages = Path(sysconfig.get_paths()["purelib"]).resolve()
    import atst_tools
    from atst_tools.calculators.abacuslite_backend import BACKEND_SOURCE, import_generalio

    atst_path = ensure_installed_origin(atst_tools.__file__, site_packages, root)
    if BACKEND_SOURCE != backend:
        raise RuntimeError(f"Expected {backend} backend, got {BACKEND_SOURCE}")
    if find_spec("abacuslite") is None:
        raise RuntimeError("ATST requires installed abacuslite")
    import abacuslite

    plugin = {
        "version": metadata.version("abacuslite"),
        "path": str(ensure_installed_origin(abacuslite.__file__, site_packages, root)),
        "files_sha256": _installed_file_digest(metadata.distribution("abacuslite"),
                                               "abacuslite/"),
        "wheel_sha256": _installed_wheel_digest(metadata.distribution("abacuslite")),
    }
    if import_generalio().__name__ != "abacuslite.io.generalio":
        raise RuntimeError("External generalio was shadowed by a legacy path")
    atst_dist = metadata.distribution("atst-tools")
    atst_wheel_sha256 = _installed_wheel_digest(atst_dist)
    if len(atst_wheel_sha256 or "") != 64:
        raise RuntimeError("ATST must be installed from a hashed wheel archive")
    distributions = {
        dist.metadata["Name"].lower().replace("_", "-"): dist.version
        for dist in metadata.distributions()
        if dist.metadata["Name"].lower().replace("_", "-") != "abacuslite"
    }
    return {
        "python": sys.version.split()[0],
        "site_packages": str(site_packages),
        "atst": {"version": metadata.version("atst-tools"), "path": str(atst_path)},
        "atst_files_sha256": _installed_file_digest(atst_dist, "atst_tools/"),
        "atst_wheel_sha256": atst_wheel_sha256,
        "distributions": dict(sorted(distributions.items())),
        "abacuslite": plugin,
    }


_GUARD_PROGRAM = """\
import builtins
from pathlib import Path
import runpy
import sys

original = builtins.__import__
def checked(name, *args, **kwargs):
    if name.startswith('atst_tools.external.ASE_interface.abacuslite'):
        raise AssertionError('Vendored backend import: ' + name)
    return original(name, *args, **kwargs)
builtins.__import__ = checked
checks = runpy.run_path(sys.argv[1])
checks['_exercise_external_consumers'](Path(sys.argv[2]))
assert not any(name.startswith('atst_tools.external.ASE_interface.abacuslite')
               for name in sys.modules)
"""


_PYTEST_PROVENANCE_PLUGIN = """\
import json
import os
from pathlib import Path
import sysconfig

def _attest(phase):
    import atst_tools
    from atst_tools.calculators.abacuslite_backend import BACKEND_SOURCE
    site = Path(sysconfig.get_paths()['purelib']).resolve()
    package = Path(atst_tools.__file__).resolve()
    root = Path(os.environ['ATST_RUNNER_ROOT']).resolve()
    expected = os.environ['ATST_RUNNER_BACKEND']
    if not package.is_relative_to(site) or package.is_relative_to(root):
        raise RuntimeError('pytest imported ATST outside site-packages: ' + str(package))
    if BACKEND_SOURCE != expected:
        raise RuntimeError(f'pytest backend {BACKEND_SOURCE} != {expected}')
    import abacuslite
    plugin = Path(abacuslite.__file__).resolve()
    if not plugin.is_relative_to(site) or plugin.is_relative_to(root):
        raise RuntimeError('pytest imported plugin outside site-packages: ' + str(plugin))
    path = Path(os.environ['ATST_RUNNER_ATTEST'])
    if phase == 'start':
        record = {'deselected': []}
    else:
        record = json.loads(path.read_text()) if path.exists() else {}
    record.setdefault('deselected', [])
    record[phase] = {'atst': str(package), 'abacuslite': str(plugin),
                     'backend': BACKEND_SOURCE}
    path.write_text(json.dumps(record, sort_keys=True))

def pytest_deselected(items):
    path = Path(os.environ['ATST_RUNNER_ATTEST'])
    record = json.loads(path.read_text()) if path.exists() else {}
    record.setdefault('deselected', []).extend(item.nodeid for item in items)
    path.write_text(json.dumps(record, sort_keys=True))

def pytest_sessionstart(session):
    _attest('start')

def pytest_sessionfinish(session, exitstatus):
    _attest('finish')
"""


_CHILD_PROVENANCE_GUARD = """\
from importlib.util import find_spec
import os
from pathlib import Path
import sysconfig

if os.environ.get('ATST_TEST_INSTALLED_WHEEL') == '1':
    site = Path(sysconfig.get_paths()['purelib']).resolve()
    root = Path(os.environ['ATST_RUNNER_ROOT']).resolve()
    package_spec = find_spec('atst_tools')
    package = Path(package_spec.origin).resolve() if package_spec else None
    if package is None:
        raise SystemExit('ATST installed package is missing')
    if not package.is_relative_to(site) or package.is_relative_to(root):
        raise SystemExit('ATST source shadow: ' + str(package))
    plugin_spec = find_spec('abacuslite')
    plugin = Path(plugin_spec.origin).resolve() if plugin_spec else None
    if plugin is None:
        raise SystemExit('abacuslite installed package is missing')
    if not plugin.is_relative_to(site) or plugin.is_relative_to(root):
        raise SystemExit('abacuslite source shadow: ' + str(plugin))
"""


def _guard_runtime_consumers(root: Path) -> None:
    """Reject a vendored import before the first ATST runtime consumer import."""
    with tempfile.TemporaryDirectory(prefix="atst-consumer-guard-") as scratch:
        result = subprocess.run(
            [sys.executable, "-c", _GUARD_PROGRAM,
             str(root / "tests/unit/test_backend_routing.py"), scratch],
            cwd=scratch, text=True, capture_output=True, timeout=90,
        )
    if result.returncode:
        raise RuntimeError("External runtime consumer imported vendored code:\n"
                           + result.stdout + result.stderr)


def _fixture_facts(root: Path) -> dict:
    """Collect numeric handoff facts checked against curated raw fixtures."""
    import numpy as np
    from atst_tools.calculators.abacuslite_backend import Abacus, import_latestio

    core = importlib.import_module(Abacus.__module__)
    fixtures = root / "src/atst_tools/external/ASE_interface/abacuslite/io/testfiles"
    with tempfile.TemporaryDirectory(prefix="atst-consumer-facts-") as scratch:
        scratch = Path(scratch)
        scf_dir = scratch / "scf"
        scf_out = scf_dir / "OUT.ABACUS"
        scf_out.mkdir(parents=True)
        scf_gold = fixtures / "multiframe_scf_trial_last"
        shutil.copy2(scf_gold / "STRU", scf_dir / "STRU")
        shutil.copy2(scf_gold / "running_scf.log", scf_out / "running_scf.log")
        template = core.AbacusTemplate()
        template.suffix = "ABACUS"
        template.calculation = "scf"
        template.atomorder = [0, 2, 1]
        scf = template.read_results(scf_dir)
        scf_forces = np.asarray(scf["forces"], dtype=float)
        assert math.isclose(float(scf["energy"]), -2898.0,
                            rel_tol=0, abs_tol=1e-10)
        np.testing.assert_allclose(
            scf_forces,
            [[-0.5, -0.01, 0], [0, 0, 0], [0.5, 0.01, 0]], atol=1e-12,
        )
        assert math.isclose(float(scf["efermi"]), -2.04,
                            rel_tol=0, abs_tol=1e-10)

        md_out = scratch / "md"
        md_out.mkdir()
        for name in ("running_md.log", "MD_dump", "eig_occ.txt"):
            shutil.copy2(fixtures / "multiframe_md_latest" / name, md_out / name)
        frames = import_latestio().read_abacus_out(md_out / "running_md.log")
        assert len(frames) == 2
        md_last = frames[-1]
        md_energy = float(md_last.get_potential_energy())
        assert math.isclose(md_energy, -1940.6311063727,
                            rel_tol=0, abs_tol=1e-9)
        np.testing.assert_allclose(
            md_last.get_forces(),
            [[0.1191299488, 0.2678988387, 1.5033416041],
             [-0.1191299488, -0.2678988387, -1.5033416041]],
            rtol=0, atol=1e-10,
        )
        np.testing.assert_allclose(md_last.cell.array, _MD_CELL_ANGSTROM,
                                   rtol=0, atol=1e-10)
        assert md_last.pbc.all()
        return {
            "scf": {
                "energy_eV": float(scf["energy"]),
                "forces_eV_per_ang": scf_forces.tolist(),
                "efermi_eV": float(scf["efermi"]),
            },
            "md": {
                "frames": len(frames),
                "last_energy_eV": md_energy,
                "last_forces_eV_per_ang": md_last.get_forces().tolist(),
                "cell_angstrom": md_last.cell.array.tolist(),
                "pbc": md_last.pbc.tolist(),
                "symbols": md_last.get_chemical_symbols(),
            },
        }


def _pytest_result(xml_path: Path) -> dict:
    tree = ET.parse(xml_path)
    suite = tree.getroot().find("testsuite")
    if suite is None:
        raise RuntimeError(f"Missing pytest testsuite in {xml_path}")
    counts = {key: int(suite.get(key, "0")) for key in
              ("tests", "failures", "errors", "skipped")}
    reasons = Counter(
        skipped.get("message", "unknown")
        for skipped in suite.iter("skipped")
    )
    testcases = {}
    for case in suite.iter("testcase"):
        node = f"{case.get('classname', '')}::{case.get('name', '')}"
        if node in testcases:
            raise RuntimeError(f"Duplicate JUnit testcase {node} in {xml_path}")
        skipped = case.find("skipped")
        if skipped is not None:
            testcases[node] = {"outcome": "skipped", "reason": skipped.get("message", "unknown")}
        elif case.find("failure") is not None:
            testcases[node] = {"outcome": "failed"}
        elif case.find("error") is not None:
            testcases[node] = {"outcome": "error"}
        else:
            testcases[node] = {"outcome": "passed"}
    if len(testcases) != counts["tests"]:
        raise RuntimeError(f"JUnit testcase count differs in {xml_path}")
    return {
        "collected": counts["tests"],
        "passed": counts["tests"] - counts["failures"] - counts["errors"] - counts["skipped"],
        "skipped": counts["skipped"],
        "failed": counts["failures"],
        "errors": counts["errors"],
        "skip_reasons": dict(sorted(reasons.items())),
        "testcases": dict(sorted(testcases.items())),
    }


def _pytest_environment(root: Path, backend: str, attest_path: Path,
                        plugin_dir: str) -> dict[str, str]:
    """Keep runner provenance while removing inherited pytest selection controls."""
    env = {**os.environ, "PYTHONPATH": plugin_dir,
           "ATST_TEST_INSTALLED_WHEEL": "1",
           "ATST_RUNNER_ROOT": str(root), "ATST_RUNNER_BACKEND": backend,
           "ATST_RUNNER_ATTEST": str(attest_path)}
    for key in ("PYTEST_ADDOPTS", "PYTEST_PLUGINS", "PYTEST_DISABLE_PLUGIN_AUTOLOAD"):
        env.pop(key, None)
    return env


def run_consumer(root: Path, backend: str, report_path: Path) -> dict:
    """Run installed-wheel provenance, guarded consumers, fixtures and tests."""
    root = root.resolve()
    if not (root / "tests").is_dir():
        raise RuntimeError(f"ATST tests directory does not exist: {root}")
    report_path = report_path.resolve()
    report_path.parent.mkdir(parents=True, exist_ok=True)
    packages = _provenance(root, backend)
    _guard_runtime_consumers(root)
    facts = _fixture_facts(root)
    xml_path = report_path.with_suffix(".pytest.xml")
    log_path = report_path.with_suffix(".pytest.log")
    attest_path = report_path.with_suffix(".pytest.attest.json")
    command = [sys.executable, "-m", "pytest", str(root / "tests"),
               "-p", "atst_installed_guard", "-o", "addopts=", "-q", "-ra",
               f"--junitxml={xml_path}"]
    with tempfile.TemporaryDirectory(prefix="atst-pytest-attest-") as plugin_dir:
        Path(plugin_dir, "atst_installed_guard.py").write_text(
            _PYTEST_PROVENANCE_PLUGIN, encoding="utf-8")
        Path(plugin_dir, "sitecustomize.py").write_text(
            _CHILD_PROVENANCE_GUARD, encoding="utf-8")
        env = _pytest_environment(root, backend, attest_path, plugin_dir)
        with log_path.open("w", encoding="utf-8") as log:
            # Some governance/example tests intentionally read repository-relative assets.
            # The temporary pytest plugin checks wheel provenance at startup and finish.
            completed = subprocess.run(command, cwd=root, env=env, stdout=log,
                                       stderr=subprocess.STDOUT, timeout=600)
    attest = json.loads(attest_path.read_text(encoding="utf-8")) if attest_path.exists() else {}
    if set(attest) != {"start", "finish", "deselected"}:
        raise RuntimeError(f"pytest wheel provenance attestation missing; see {log_path}")
    tests = _pytest_result(xml_path)
    tests["deselected"] = len(attest["deselected"])
    report = {
        "schema": SCHEMA,
        "backend": backend,
        "packages": packages,
        "tests": tests,
        "facts": facts,
        "pytest_log": str(log_path),
        "pytest_xml": str(xml_path),
        "pytest_attestation": attest,
        "pytest_exit_code": completed.returncode,
        "pytest_cwd": str(root),
        "pytest_command": command,
    }
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    if completed.returncode or tests["failed"] or tests["errors"]:
        raise RuntimeError(f"ATST installed consumer suite failed; see {log_path}")
    verify_report(report)
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atst-root", type=Path)
    parser.add_argument("--backend", choices=("external",))
    parser.add_argument("--report", type=Path)
    parser.add_argument("--verify-report", type=Path)
    args = parser.parse_args(argv)
    if args.verify_report:
        with args.verify_report.open(encoding="utf-8") as stream:
            verify_report(json.load(stream))
        print("ATST external consumer report verified")
        return 0
    if not all((args.atst_root, args.backend, args.report)):
        parser.error("--atst-root, --backend and --report are required together")
    report = run_consumer(args.atst_root, args.backend, args.report)
    summary = {key: value for key, value in report["tests"].items()
               if key != "testcases"}
    print(json.dumps({"backend": args.backend, "tests": summary,
                      "report": str(args.report)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
