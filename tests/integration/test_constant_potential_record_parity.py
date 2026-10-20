"""Opt-in parity check for FCP v2 and ATST constant-potential records."""

from __future__ import annotations

import importlib.util
import os
from pathlib import Path
import subprocess
import sys

import pytest


pytestmark = [pytest.mark.integration, pytest.mark.slurm]

ROOT = Path(__file__).resolve().parents[2]
PRODUCER = ROOT / "scripts" / "run_fcp_v2_reference.py"
EXTRACTOR = ROOT / "scripts" / "extract_constant_potential_record.py"
COMPARATOR = ROOT / "scripts" / "compare_constant_potential_results.py"


def _required_env(name: str) -> str:
    value = os.environ.get(name)
    if not value:
        pytest.skip(f"set {name} for constant-potential record parity")
    return value


def test_fcp_v2_record_matches_atst_record(tmp_path):
    """Compare one real FCP v2 record with one real ATST record.

    Required environment:

    - ``ATST_FCP_REFERENCE_CONFIG``: producer JSON configuration.
    - ``ATST_FCP_REFERENCE_SOURCE``: pinned ``version2/FCPelectrochem.py``.
    - ``ATST_CP_CANDIDATE_RESULTS``: ATST ``constant_potential_results.json``.
    - ``ATST_CP_CANDIDATE_POINT``: optional zero-based point index.
    """
    config = _required_env("ATST_FCP_REFERENCE_CONFIG")
    source = _required_env("ATST_FCP_REFERENCE_SOURCE")
    candidate_results = _required_env("ATST_CP_CANDIDATE_RESULTS")
    point = os.environ.get("ATST_CP_CANDIDATE_POINT", "0")
    if importlib.util.find_spec("ase.calculators.abacus") is None:
        pytest.skip("legacy ase-abacus is not installed")

    reference_record = tmp_path / "fcp_reference.json"
    candidate_record = tmp_path / "atst_candidate.json"

    produced = subprocess.run(
        [
            sys.executable,
            str(PRODUCER),
            "--config",
            config,
            "--fcp-source",
            source,
            "--output",
            str(reference_record),
        ],
        capture_output=True,
        text=True,
        timeout=7200,
    )
    assert produced.returncode == 0, produced.stdout + produced.stderr

    extracted = subprocess.run(
        [
            sys.executable,
            str(EXTRACTOR),
            "--input",
            candidate_results,
            "--output",
            str(candidate_record),
            "--point",
            point,
        ],
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert extracted.returncode == 0, extracted.stdout + extracted.stderr

    compared = subprocess.run(
        [
            sys.executable,
            str(COMPARATOR),
            "--reference",
            str(reference_record),
            "--candidate",
            str(candidate_record),
        ],
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert compared.returncode == 0, compared.stdout + compared.stderr
