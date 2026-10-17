"""Subprocesses in installed-wheel parity tests must use the same wheel."""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import sysconfig

import pytest

from helpers import INSTALLED_WHEEL_TEST_ENV, child_test_environment


SOURCE_ROOT = Path(__file__).resolve().parents[2] / "src"


def _child_package_path(environment: dict[str, str], cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-c", "import atst_tools; print(atst_tools.__file__)"],
        cwd=cwd, env=environment, capture_output=True, text=True, timeout=30,
    )


def test_child_environment_preserves_installed_guard_and_rejects_source_path():
    guarded = {INSTALLED_WHEEL_TEST_ENV: "1", "PYTHONPATH": "/temporary/attest",
               "CUSTOM_RUNTIME_DEP": "kept"}
    assert child_test_environment(SOURCE_ROOT, base=guarded) == guarded
    shadowed = {**guarded, "PYTHONPATH": str(SOURCE_ROOT)}
    with pytest.raises(AssertionError, match="source path"):
        child_test_environment(SOURCE_ROOT, base=shadowed)


def test_child_import_uses_selected_source_and_rejects_installed_shadow(tmp_path):
    environment = child_test_environment(SOURCE_ROOT)
    actual = _child_package_path(environment, tmp_path)
    assert actual.returncode == 0, actual.stderr
    package = Path(actual.stdout.strip()).resolve()
    if os.environ.get(INSTALLED_WHEEL_TEST_ENV) == "1":
        assert package.is_relative_to(Path(sysconfig.get_paths()["purelib"]).resolve())
        shadowed = dict(environment)
        shadowed["PYTHONPATH"] = str(SOURCE_ROOT) + os.pathsep + environment.get("PYTHONPATH", "")
        rejected = _child_package_path(shadowed, tmp_path)
        assert rejected.returncode != 0
        assert "ATST source shadow" in rejected.stderr
    else:
        assert package.is_relative_to(SOURCE_ROOT)
