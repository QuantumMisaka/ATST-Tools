"""Consumer regressions for a single ABACUS backend decision."""

from __future__ import annotations

import builtins
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

import pytest


def _reject_vendored_imports(monkeypatch):
    original_import = builtins.__import__

    def checked_import(name, *args, **kwargs):
        if name.startswith("atst_tools.external.ASE_interface.abacuslite"):
            raise AssertionError(f"Vendored backend import: {name}")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", checked_import)


def test_external_consumers_do_not_import_vendored_backend(tmp_path):
    """Guard the first consumer import, independent of pytest collection order."""
    script = """
import runpy
import sys
from pathlib import Path
import pytest

checks = runpy.run_path(sys.argv[1])
with pytest.MonkeyPatch.context() as patch:
    checks['_reject_vendored_imports'](patch)
    checks['_exercise_external_consumers'](Path(sys.argv[2]))
assert not any(name.startswith('atst_tools.external.ASE_interface.abacuslite')
               for name in sys.modules)
"""
    result = subprocess.run(
        [sys.executable, "-c", script, str(Path(__file__).resolve()), str(tmp_path)],
        capture_output=True, text=True, timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def _exercise_external_consumers(tmp_path):
    """Exercise actual structure, reverse, MD and gate consumer entry points."""
    from atst_tools.calculators.abacuslite_backend import BACKEND_SOURCE

    assert BACKEND_SOURCE == "external"

    from atst_tools.utils import abacus_io, gate_compensation, io, reverse_config
    from atst_tools.workflows import md

    assert abacus_io._import_generalio().__name__ == "abacuslite.io.generalio"
    assert abacus_io._import_latestio().__module__ == "abacuslite.io.latestio"
    assert md._import_legacyio().__module__ == "abacuslite.io.legacyio"
    assert reverse_config.endpoint_has_energy_forces(tmp_path) is False

    run = tmp_path / "run"
    run.mkdir()
    (run / "INPUT").write_text("INPUT_PARAMETERS\ngamma_only 1\n", encoding="utf-8")
    (run / "STRU").write_text(
        "ATOMIC_SPECIES\nH 1 H.upf\nLATTICE_CONSTANT\n1\n"
        "LATTICE_VECTORS\n10 0 0\n0 10 0\n0 0 10\n"
        "ATOMIC_POSITIONS\nDirect\nH\n0\n1\n.2 .2 .2 1 1 1\n",
        encoding="utf-8",
    )
    assert io.read_abacus_stru(run / "STRU").get_chemical_symbols() == ["H"]
    assert reverse_config._read_input(run)["gamma_only"] == "1"
    assert reverse_config._read_stru_species(run)[0] == {"H": "H.upf"}
    assert reverse_config._resolve_kpts(run, {"gamma_only": "1"}) == [1, 1, 1]
    assert reverse_config._derive_kgrid_from_kspacing(run, "1") is not None

    # A missing cube is an expected file error; the parser import must already
    # have been resolved from the selected external backend.
    with pytest.raises((ValueError, FileNotFoundError)):
        gate_compensation.read_gate_compensation(run, 1.0, 1.0)


def test_missing_internal_external_module_never_uses_vendored(monkeypatch):
    from atst_tools.calculators import abacuslite_backend

    assert abacuslite_backend.BACKEND_SOURCE == "external"
    original_import = builtins.__import__

    def broken_generalio(name, *args, **kwargs):
        if name == "abacuslite.io":
            raise ModuleNotFoundError("broken installed plugin", name="abacuslite.io")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", broken_generalio)
    from atst_tools.utils import abacus_io

    with pytest.raises(ModuleNotFoundError, match="broken installed plugin"):
        abacus_io._import_generalio()


def test_external_frame_selection_fermi_order_and_native_md(monkeypatch, tmp_path):
    """Use real parser fixtures through installed calculator and MD consumers."""
    from atst_tools.calculators.abacuslite_backend import BACKEND_SOURCE

    assert BACKEND_SOURCE == "external"
    from atst_tools.workflows.md import AbacusNativeMDRunner

    fixtures = (
        Path(__file__).resolve().parents[2]
        / "src/atst_tools/external/ASE_interface/abacuslite/io/testfiles"
    )
    gold = fixtures / "multiframe_scf_trial_last"
    run = tmp_path / "scf"
    out = run / "OUT.ABACUS"
    out.mkdir(parents=True)
    shutil.copy2(gold / "STRU", run / "STRU")
    shutil.copy2(gold / "running_scf.log", out / "running_scf.log")

    _reject_vendored_imports(monkeypatch)
    from abacuslite.core import AbacusTemplate

    template = AbacusTemplate()
    template.suffix = "ABACUS"
    template.calculation = "scf"
    template.atomorder = [0, 2, 1]
    results = template.read_results(run)
    np.testing.assert_allclose(
        results["forces"],
        [[-0.5, -0.01, 0], [0, 0, 0], [0.5, 0.01, 0]],
        atol=1e-12,
    )
    assert results["efermi"] == pytest.approx(-2.04)
    assert results["fermi_level"] == pytest.approx(results["efermi"])

    md_run = tmp_path / "md"
    md_output = md_run / "OUT.ABACUS"
    md_output.mkdir(parents=True)
    (md_run / "INPUT").write_text("INPUT_PARAMETERS\nsuffix ABACUS\n", encoding="utf-8")
    for name in ("running_md.log", "MD_dump", "eig_occ.txt"):
        shutil.copy2(fixtures / "multiframe_md_latest" / name, md_output / name)
    runner = AbacusNativeMDRunner({}, "abacus", {"directory": str(md_run)})
    frames = runner._collect_frames()
    assert len(frames) == 2
    assert frames[-1].get_chemical_symbols() == ["As", "Ga"]
    assert frames[-1].pbc.all()
    np.testing.assert_allclose(frames[-1].cell.array[0, 0], 5.34197 * 0.529177, atol=1e-10)
    assert np.isfinite(frames[-1].get_potential_energy())

    from atst_tools.utils.abacus_io import collect_abacus_output

    summary = collect_abacus_output(str(md_run), str(tmp_path / "summary.json"))
    assert summary["backend"] == "external"
    assert summary["parsed"] is True
    assert summary["energy_eV"] == pytest.approx(-1940.6311063727)
    assert summary["max_force_eV_per_ang"] > 0


def test_missing_plugin_fails_closed_before_a_vendored_import():
    """A missing required package cannot revive the historical snapshot."""
    script = """
import builtins
from importlib import util
original_find_spec = util.find_spec
original_import = builtins.__import__
def missing_plugin(name, *args, **kwargs):
    if name == 'abacuslite':
        return None
    return original_find_spec(name, *args, **kwargs)
def reject_vendored(name, *args, **kwargs):
    if name.startswith('atst_tools.external.ASE_interface.abacuslite'):
        raise AssertionError('legacy implementation imported')
    return original_import(name, *args, **kwargs)
util.find_spec = missing_plugin
builtins.__import__ = reject_vendored
from atst_tools.calculators import abacuslite_backend
abacuslite_backend._load_abacuslite_backend()
"""
    result = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True)
    assert result.returncode != 0
    assert "abacuslite is required" in result.stderr
    assert "legacy implementation imported" not in result.stderr


def test_legacy_import_paths_forward_to_the_installed_modules():
    """External callers retain imports without a second parser or calculator."""
    import abacuslite.core as core
    from abacuslite.io import generalio, latestio, legacyio
    from abacuslite.utils import ksampling
    from atst_tools.external.ASE_interface.abacuslite import Abacus
    from atst_tools.external.ASE_interface.abacuslite import core as old_core
    from atst_tools.external.ASE_interface.abacuslite.io import (
        generalio as old_generalio, latestio as old_latestio, legacyio as old_legacyio,
    )
    from atst_tools.external.ASE_interface.abacuslite.utils import ksampling as old_ksampling

    assert Abacus is core.Abacus
    assert (old_core, old_generalio, old_latestio, old_legacyio, old_ksampling) == (
        core, generalio, latestio, legacyio, ksampling,
    )
