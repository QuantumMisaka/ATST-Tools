"""Reproduce the committed FCP v2 trace from its pinned source."""

from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest


pytestmark = pytest.mark.integration

ROOT = Path(__file__).resolve().parents[2]
REFERENCE = ROOT / "tests" / "data" / "constant_potential_reference.json"
GENERATOR = ROOT / "tests" / "data" / "generate_constant_potential_reference.py"


def test_committed_fcp_v2_reference_is_reproducible(tmp_path):
    """Regenerate the reference trace without touching the checked-out tree.

    Set ``ATST_FCP_REFERENCE_SOURCE`` to the pinned
    ``version2/FCPelectrochem.py`` checkout. The generator verifies its SHA256
    and uses a deterministic analytic backend; it does not execute ABACUS.
    """
    source = os.environ.get("ATST_FCP_REFERENCE_SOURCE")
    if not source:
        pytest.skip(
            "set ATST_FCP_REFERENCE_SOURCE to the pinned v2_FCPelectrochem.py"
        )
    pytest.importorskip("pandas")

    source_path = Path(source).expanduser().resolve()
    if not source_path.is_file():
        pytest.fail(f"ATST_FCP_REFERENCE_SOURCE is not a file: {source_path}")

    shutil.copy2(GENERATOR, tmp_path / "generate.py")
    shutil.copy2(REFERENCE, tmp_path / "constant_potential_reference.json")
    result = subprocess.run(
        [sys.executable, str(tmp_path / "generate.py"), str(source_path)],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    actual = json.loads(
        (tmp_path / "constant_potential_reference.json").read_text(encoding="utf-8")
    )
    expected = json.loads(REFERENCE.read_text(encoding="utf-8"))
    assert actual == expected
