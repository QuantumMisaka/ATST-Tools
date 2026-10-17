#!/usr/bin/env python3
"""Run the same installed-wheel ATST consumer suite for either abacuslite backend.

The caller installs wheels in separate environments, pins their source commits,
and invokes this script from outside both checkouts. No engine or network is used.
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


SCHEMA = "atst-abacuslite-consumer-v1"
_BACKEND_SKIP_FRAGMENTS = (
    "external-package lane",
    "no-plugin environment",
)


def ensure_installed_origin(path: str | Path, site_packages: Path, checkout: Path) -> Path:
    """Require a wheel import from this interpreter's site-packages directory."""
    resolved = Path(path).resolve()
    site_packages = site_packages.resolve()
    checkout = checkout.resolve()
    if resolved.is_relative_to(checkout) or not resolved.is_relative_to(site_packages):
        raise RuntimeError(f"Expected installed site-packages import, got {resolved}")
    return resolved


def _non_backend_skips(reasons: dict[str, int]) -> dict[str, int]:
    return {
        reason: count for reason, count in reasons.items()
        if not any(fragment in reason for fragment in _BACKEND_SKIP_FRAGMENTS)
    }


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


def compare_reports(vendored: dict, external: dict) -> dict:
    """Reject mismatched collection, unexpected skips, or fixture facts."""
    for expected, report in (("vendored", vendored), ("external", external)):
        if report.get("schema") != SCHEMA or report.get("backend") != expected:
            raise ValueError(f"Expected {expected} {SCHEMA} report")
        tests = report.get("tests", {})
        if tests.get("collected", 0) <= 0:
            raise ValueError(f"{expected} collected no tests")
        packages = report.get("packages", {})
        if not packages.get("python") or not packages.get("distributions"):
            raise ValueError(f"{expected} lacks installed dependency provenance")
        for name in ("atst_files_sha256", "atst_wheel_sha256"):
            if len(packages.get(name) or "") != 64:
                raise ValueError(f"{expected} lacks {name} provenance")
        if tests.get("failed") or tests.get("errors"):
            raise ValueError(f"{expected} offline suite failed")
        if tests.get("passed", 0) + tests.get("skipped", 0) != tests.get("collected"):
            raise ValueError(f"{expected} test collection is incomplete")
    if vendored["packages"]["atst"]["version"] != external["packages"]["atst"]["version"]:
        raise ValueError("ATST wheel versions differ across backends")
    for name in ("python", "atst_files_sha256", "atst_wheel_sha256", "distributions"):
        if vendored["packages"].get(name) != external["packages"].get(name):
            raise ValueError(f"Installed {name} differs across backends")
    if vendored["tests"]["collected"] != external["tests"]["collected"]:
        raise ValueError("ATST test collection differs across backends")
    vendored_cases = vendored["tests"].get("testcases", {})
    external_cases = external["tests"].get("testcases", {})
    if not vendored_cases or vendored_cases.keys() != external_cases.keys():
        raise ValueError("ATST test case identities differ across backends")
    for node, left in vendored_cases.items():
        right = external_cases[node]
        if left == right:
            continue
        if any(case.get("outcome") == "skipped" and
               any(fragment in case.get("reason", "") for fragment in _BACKEND_SKIP_FRAGMENTS)
               for case in (left, right)):
            continue
        raise ValueError(f"ATST test outcome differs for {node}: {left} != {right}")
    if _non_backend_skips(vendored["tests"]["skip_reasons"]) != _non_backend_skips(
        external["tests"]["skip_reasons"]
    ):
        raise ValueError("Non-backend ATST skip reasons differ across backends")
    _compare_values(vendored["facts"], external["facts"], "facts")
    return {}


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
    if Path.cwd().resolve().is_relative_to(root):
        raise RuntimeError("Run installed consumer checks from outside the ATST checkout")
    site_packages = Path(sysconfig.get_paths()["purelib"]).resolve()
    import atst_tools
    from atst_tools.calculators.abacuslite_backend import BACKEND_SOURCE, import_generalio

    atst_path = ensure_installed_origin(atst_tools.__file__, site_packages, root)
    if BACKEND_SOURCE != backend:
        raise RuntimeError(f"Expected {backend} backend, got {BACKEND_SOURCE}")
    plugin = None
    if backend == "external":
        if find_spec("abacuslite") is None:
            raise RuntimeError("External backend requested without installed abacuslite")
        import abacuslite

        plugin = {
            "version": metadata.version("abacuslite"),
            "path": str(ensure_installed_origin(abacuslite.__file__, site_packages, root)),
            "files_sha256": _installed_file_digest(metadata.distribution("abacuslite"),
                                                   "abacuslite/"),
            "wheel_sha256": _installed_wheel_digest(metadata.distribution("abacuslite")),
        }
        if import_generalio().__name__ != "abacuslite.io.generalio":
            raise RuntimeError("External generalio was shadowed by the vendored snapshot")
    elif find_spec("abacuslite") is not None:
        raise RuntimeError("Vendored backend requires an environment without abacuslite")
    elif not import_generalio().__name__.startswith("atst_tools.external.ASE_interface.abacuslite"):
        raise RuntimeError("Vendored generalio did not originate from ATST")
    atst_dist = metadata.distribution("atst-tools")
    atst_wheel_sha256 = _installed_wheel_digest(atst_dist)
    if len(atst_wheel_sha256 or "") != 64:
        raise RuntimeError("ATST must be installed from a hashed wheel archive")
    if plugin is not None and len(plugin["wheel_sha256"] or "") != 64:
        raise RuntimeError("abacuslite must be installed from a hashed wheel archive")
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
    if expected == 'external':
        import abacuslite
        plugin = Path(abacuslite.__file__).resolve()
        if not plugin.is_relative_to(site) or plugin.is_relative_to(root):
            raise RuntimeError('pytest imported plugin outside site-packages: ' + str(plugin))
    else:
        plugin = None
    path = Path(os.environ['ATST_RUNNER_ATTEST'])
    record = json.loads(path.read_text()) if path.exists() else {}
    record[phase] = {'atst': str(package), 'abacuslite': str(plugin) if plugin else None,
                     'backend': BACKEND_SOURCE}
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
    if os.environ['ATST_RUNNER_BACKEND'] == 'external':
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
        np.testing.assert_allclose(md_last.cell.array[0, 0], 5.34197 * 0.529177,
                                   atol=1e-10)
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


def run_consumer(root: Path, backend: str, report_path: Path) -> dict:
    """Run installed-wheel provenance, guarded consumers, fixtures and tests."""
    root = root.resolve()
    if not (root / "tests").is_dir():
        raise RuntimeError(f"ATST tests directory does not exist: {root}")
    report_path = report_path.resolve()
    report_path.parent.mkdir(parents=True, exist_ok=True)
    packages = _provenance(root, backend)
    if backend == "external":
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
        env = {**os.environ, "PYTHONPATH": plugin_dir,
               "ATST_TEST_INSTALLED_WHEEL": "1",
               "ATST_RUNNER_ROOT": str(root), "ATST_RUNNER_BACKEND": backend,
               "ATST_RUNNER_ATTEST": str(attest_path)}
        with log_path.open("w", encoding="utf-8") as log:
            # Some governance/example tests intentionally read repository-relative assets.
            # The temporary pytest plugin checks wheel provenance at startup and finish.
            completed = subprocess.run(command, cwd=root, env=env, stdout=log,
                                       stderr=subprocess.STDOUT, timeout=600)
    attest = json.loads(attest_path.read_text(encoding="utf-8")) if attest_path.exists() else {}
    if set(attest) != {"start", "finish"}:
        raise RuntimeError(f"pytest wheel provenance attestation missing; see {log_path}")
    tests = _pytest_result(xml_path)
    report = {
        "schema": SCHEMA,
        "backend": backend,
        "packages": packages,
        "tests": tests,
        "facts": facts,
        "pytest_log": str(log_path),
        "pytest_xml": str(xml_path),
        "pytest_attestation": attest,
        "pytest_cwd": str(root),
        "pytest_command": command,
    }
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    if completed.returncode or tests["failed"] or tests["errors"]:
        raise RuntimeError(f"ATST installed consumer suite failed; see {log_path}")
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atst-root", type=Path)
    parser.add_argument("--backend", choices=("vendored", "external"))
    parser.add_argument("--report", type=Path)
    parser.add_argument("--compare", nargs=2, metavar=("VENDORED_JSON", "EXTERNAL_JSON"),
                        type=Path)
    args = parser.parse_args(argv)
    if args.compare:
        with args.compare[0].open(encoding="utf-8") as stream:
            vendored = json.load(stream)
        with args.compare[1].open(encoding="utf-8") as stream:
            external = json.load(stream)
        compare_reports(vendored, external)
        print("ATST vendored/external consumer parity passed")
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
