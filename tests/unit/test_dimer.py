"""Dimer optimizer facts and artifact handoff tests."""

import json

import pytest
from ase import Atoms

from atst_tools.mep import dimer as dimer_module


@pytest.mark.parametrize("signal, expected", [(True, True), (False, False), (None, None), (0, None)])
def test_dimer_persists_optimizer_owned_convergence(monkeypatch, tmp_path, capsys, signal, expected):
    """Completed optimization remains complete even when unconverged."""

    class FakeControl:
        def __init__(self, **kwargs):
            pass

    class FakeMinModeAtoms:
        def __init__(self, atoms, control):
            pass

        def displace(self, displacement_vector=None):
            pass

    class FakeTranslate:
        def __init__(self, atoms):
            self.nsteps = 7

        def attach(self, observer):
            pass

        def run(self, fmax=None, steps=None):
            return signal

    class FakeTrajectory:
        def __init__(self, *args):
            pass

        def close(self):
            pass

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(dimer_module, "DimerControl", FakeControl)
    monkeypatch.setattr(dimer_module, "MinModeAtoms", FakeMinModeAtoms)
    monkeypatch.setattr(dimer_module, "MinModeTranslate", FakeTranslate)
    monkeypatch.setattr(dimer_module, "Trajectory", FakeTrajectory)
    monkeypatch.setattr(dimer_module.AbacusDimer, "set_calculator", lambda self: None)
    dimer = dimer_module.AbacusDimer(
        Atoms("H"), {}, "abacus", {}, displacement_vector=[[0.01, 0, 0]],
    )

    dimer.run(fmax=0.07, max_steps=11)

    stage = json.loads((tmp_path / "atst_artifacts.json").read_text())["stages"][0]
    assert stage == {
        "name": "dimer", "status": "complete", "converged": expected,
        "role": "final", "criterion": "dimer", "fmax": 0.07,
        "fmax_unit": "eV/Angstrom", "steps": 11, "actual_steps": 7,
    }
    assert dimer.last_stage_record.to_manifest() == stage
    assert ("workflow=dimer" in capsys.readouterr().out) is (expected is False)
    assert not (tmp_path / "ts_validation.json").exists()
