"""Regression checks for the generated NEB band and restart contract."""

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.emt import EMT
from ase.io import read, write

from atst_tools.workflows import execution


def _run_neb_with_generated_chain(monkeypatch, tmp_path, *, ts_guess):
    """Run the public NEB executor with real path generation and ASE IO."""
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        execution.CalculatorFactory,
        "get_calculator",
        lambda *args, **kwargs: EMT(),
    )
    for name, x in (("start", 0.0), ("ts", 0.3), ("end", 0.6)):
        write(f"{name}.traj", Atoms("H", positions=[[x, 0.0, 0.0]]))
    make = {
        "n_images": 3,
        "init_structure": "start.traj",
        "final_structure": "end.traj",
        "method": "linear",
        "no_align": True,
    }
    if ts_guess:
        make["ts_guess"] = "ts.traj"
    config = {
        "type": "neb",
        "trajectory": "neb.traj",
        "parallel": False,
        "fmax": 1.0,
        "max_steps": 0,
        "make": make,
    }
    generated = execution.run_neb({}, "abacus", config)
    return config, generated


@pytest.mark.parametrize("bands", [1, 5])
def test_ts_guess_restart_restores_complete_last_band(monkeypatch, tmp_path, bands):
    config, generated = _run_neb_with_generated_chain(
        monkeypatch, tmp_path, ts_guess=True
    )
    assert len(generated) == 6
    recorded = read("neb.traj", index=":")
    assert len(recorded) == 6
    if bands > 1:
        write("neb.traj", recorded * bands)
    config["restart"] = True

    restarted = execution.run_neb({}, "abacus", config)

    assert len(restarted) == 6
    np.testing.assert_allclose(
        [image.positions[0, 0] for image in restarted],
        [image.positions[0, 0] for image in recorded],
    )
    assert restarted[0].positions[0, 0] == pytest.approx(0.0)
    assert restarted[-1].positions[0, 0] == pytest.approx(0.6)


def test_ordinary_generated_band_restores_five_images(monkeypatch, tmp_path):
    config, generated = _run_neb_with_generated_chain(
        monkeypatch, tmp_path, ts_guess=False
    )
    assert len(generated) == 5
    recorded = read("neb.traj", index=":")
    assert len(recorded) == 5
    config["restart"] = True

    restarted = execution.run_neb({}, "abacus", config)

    assert len(restarted) == 5
    np.testing.assert_allclose(
        [image.positions[0, 0] for image in restarted],
        [image.positions[0, 0] for image in recorded],
    )


@pytest.mark.parametrize("frame_count", [5, 29])
def test_ts_guess_restart_rejects_incomplete_band(
    monkeypatch, tmp_path, frame_count
):
    config, _ = _run_neb_with_generated_chain(monkeypatch, tmp_path, ts_guess=True)
    recorded = read("neb.traj", index=":")
    write("neb.traj", (recorded * 5)[:frame_count])
    config["restart"] = True

    with pytest.raises(ValueError, match="(fewer than expected|whole number of bands)"):
        execution.run_neb({}, "abacus", config)
