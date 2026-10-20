import json
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read, write

from atst_tools.calculators.constant_potential import (
    CP_FACTS_INFO_KEY,
    CP_IDENTITY_INFO_KEY,
    publish_constant_potential_facts,
)


def _atoms(energy=None, force=0.0, momenta=False):
    atoms = Atoms("H2", positions=[[0, 0, 0], [0.75, 0, 0]], cell=[8, 8, 8], pbc=True)
    if momenta:
        atoms.set_momenta([[0.1, 0.0, 0.0], [-0.1, 0.0, 0.0]])
    if energy is not None:
        atoms.calc = SinglePointCalculator(
            atoms,
            energy=energy,
            forces=np.full((len(atoms), 3), force),
        )
    return atoms


def _cp_atoms():
    """Return a saved MD frame with a validated compensated-gate fact envelope."""
    atoms = _atoms(4.2, 0.1, momenta=True)
    identity = {
        "energy_boundary": "compensated_gate",
        "potential_v": None,
        "target_mu_ev": -15.0,
        "reference_electrode": "custom",
        "work_ref": None,
        "work_ref_source": None,
        "reference_pH": None,
        "temperature_K": None,
        "reference_electrons": 10.0,
        "boundary_parameters": {},
        "energy_definition": "omega_compensated_gate",
        "calculator": "constant_potential",
    }
    publish_constant_potential_facts(
        atoms,
        {
            "energy": 4.2,
            "free_energy": 4.1,
            "raw_energy": 1.2,
            "raw_free_energy": 1.1,
            "efermi": -15.78,
            "vacuum_level": None,
            "fermishift": None,
            "energy_boundary": "compensated_gate",
            "boundary_parameters": {},
            "compensation": {
                "boundary": "compensated_gate",
                "gate_derivative_ev": 0.4,
                "dipole_derivative_ev": 0.4,
                "compensation_derivative_ev": 0.8,
                "integrated_electrons": 10.2,
                "density_electron_error": 0.0,
                "ionic_valence_electrons": 10.0,
                "total_dipole_ry_au": 0.0,
                "density_precision": 12,
                "density_electron_tolerance": 1e-8,
            },
            "mu_calc": -14.98,
            "mu_target": -15.0,
            "potential_calc": None,
            "potential_target": None,
            "residual_mu": 0.02,
            "residual_v": 0.02,
            "omega_correction": 3.0,
            "nelec": 10.2,
            "reference_electrons": 10.0,
            "delta_nelec": 0.2,
            "forces": [[0.1, 0.1, 0.1], [0.1, 0.1, 0.1]],
            "scf_converged": True,
            "cp_converged": True,
            "cp_iterations": 2,
            "cp_history": [{"converged": True}],
        },
        identity,
    )
    atoms.info["atst_md_step"] = 7
    atoms.info["atst_md_time_fs"] = 3.5
    atoms.calc = None
    return atoms


def test_summarize_md_trajectory_reports_energy_temperature_force_and_tail(tmp_path):
    from atst_tools.utils.md_post import summarize_md_trajectory

    traj = tmp_path / "md.traj"
    write(traj, [_atoms(-1.0, 0.1, momenta=True), _atoms(-1.5, 0.2, momenta=True)])

    summary = summarize_md_trajectory(traj, tail=1, metadata={"driver": "ase"})

    assert summary["workflow"] == "md"
    assert summary["status"]["n_frames"] == 2
    assert summary["latest"]["step"] == 1
    assert "time_fs" not in summary["latest"]
    assert summary["latest"]["energy_eV"] == pytest.approx(-1.5)
    assert summary["latest"]["max_force_eV_per_A"] == pytest.approx(np.linalg.norm([0.2, 0.2, 0.2]))
    assert summary["latest"]["temperature_K"] > 0
    assert summary["metadata"]["driver"] == "ase"
    assert [frame["step"] for frame in summary["frames"]] == [1]


def test_summarize_md_trajectory_uses_saved_step_time_and_valid_cp_facts(tmp_path):
    from atst_tools.utils.md_post import summarize_md_trajectory

    traj = tmp_path / "cp_md.traj"
    write(traj, _cp_atoms())

    summary = summarize_md_trajectory(traj)
    latest = summary["latest"]
    assert latest["step"] == 7
    assert latest["time_fs"] == pytest.approx(3.5)
    assert latest["constant_potential"] == {
        "energy_boundary": "compensated_gate",
        "energy_definition": "omega_compensated_gate",
        "nelec": pytest.approx(10.2),
        "mu_target_ev": pytest.approx(-15.0),
        "mu_calc_ev": pytest.approx(-14.98),
        "residual_mu_ev": pytest.approx(0.02),
        "cp_iterations": 2,
        "raw_energy_eV": pytest.approx(1.2),
        "raw_free_energy_eV": pytest.approx(1.1),
        "grand_potential_eV": pytest.approx(4.2),
        "grand_potential_free_energy_eV": pytest.approx(4.1),
    }


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("atst_md_step", -1),
        ("atst_md_step", 1.5),
        ("atst_md_step", True),
        ("atst_md_time_fs", -0.1),
        ("atst_md_time_fs", float("nan")),
        ("atst_md_time_fs", "3.5"),
        ("atst_md_time_fs", True),
    ],
)
def test_summarize_md_trajectory_rejects_invalid_saved_step_or_time(tmp_path, field, value):
    from atst_tools.utils.md_post import summarize_md_trajectory

    atoms = _atoms(-1.0)
    atoms.info.update({"atst_md_step": 0, "atst_md_time_fs": 0.0})
    atoms.info[field] = value
    traj = tmp_path / "bad_md.traj"
    write(traj, atoms)

    with pytest.raises(ValueError, match=field):
        summarize_md_trajectory(traj)


@pytest.mark.parametrize("missing", ["atst_md_step", "atst_md_time_fs"])
def test_summarize_md_trajectory_rejects_incomplete_saved_step_time(tmp_path, missing):
    from atst_tools.utils.md_post import summarize_md_trajectory

    atoms = _atoms(-1.0)
    atoms.info.update({"atst_md_step": 0, "atst_md_time_fs": 0.0})
    del atoms.info[missing]
    traj = tmp_path / "incomplete_md.traj"
    write(traj, atoms)

    with pytest.raises(ValueError, match=missing):
        summarize_md_trajectory(traj)


@pytest.mark.parametrize("damage", ["missing_mu", "stale_geometry", "orphan_identity"])
def test_summarize_md_trajectory_rejects_invalid_cp_payload(tmp_path, damage):
    from atst_tools.utils.md_post import summarize_md_trajectory

    atoms = _cp_atoms()
    if damage == "missing_mu":
        atoms.info[CP_FACTS_INFO_KEY]["facts"].pop("mu_calc")
    elif damage == "stale_geometry":
        atoms.positions[0, 0] += 0.1
    else:
        del atoms.info[CP_FACTS_INFO_KEY]
        assert CP_IDENTITY_INFO_KEY in atoms.info
    traj = tmp_path / "invalid_cp_md.traj"
    write(traj, atoms)

    with pytest.raises(ValueError, match="constant-potential facts"):
        summarize_md_trajectory(traj)


def test_summarize_md_trajectory_handles_missing_results(tmp_path):
    from atst_tools.utils.md_post import summarize_md_trajectory

    traj = tmp_path / "md.traj"
    write(traj, [_atoms()])

    summary = summarize_md_trajectory(traj)
    latest = summary["latest"]

    assert np.isnan(latest["energy_eV"])
    assert np.isnan(latest["max_force_eV_per_A"])
    assert np.isnan(latest["temperature_K"])


def test_post_md_trajectory_writes_summary_and_default_extxyz(tmp_path):
    from atst_tools.utils.md_post import post_md_trajectory

    traj = tmp_path / "md.traj"
    write(traj, [_atoms(-1.0, 0.1), _atoms(-1.5, 0.2)])

    result = post_md_trajectory(
        traj,
        output_prefix=tmp_path / "md_post",
        summary_output=tmp_path / "md_post_summary.json",
    )

    assert (tmp_path / "md_post.extxyz").exists()
    assert len(read(tmp_path / "md_post.extxyz", index=":")) == 2
    summary = json.loads((tmp_path / "md_post_summary.json").read_text(encoding="utf-8"))
    assert summary["workflow"] == "md"
    assert result["converted"]["path"] == str(tmp_path / "md_post.extxyz")


def test_post_md_trajectory_supports_frame_and_stride(tmp_path):
    from atst_tools.utils.md_post import post_md_trajectory

    traj = tmp_path / "md.traj"
    write(traj, [_atoms(-1.0), _atoms(-1.1), _atoms(-1.2), _atoms(-1.3)])

    frame_result = post_md_trajectory(
        traj,
        output_prefix=tmp_path / "last",
        output_format="xyz",
        summary_output=tmp_path / "last_summary.json",
        frame=-1,
    )
    stride_result = post_md_trajectory(
        traj,
        output_prefix=tmp_path / "stride",
        output_format="extxyz",
        summary_output=tmp_path / "stride_summary.json",
        stride=2,
    )

    assert len(read(frame_result["converted"]["path"], index=":")) == 1
    assert len(read(stride_result["converted"]["path"], index=":")) == 2
    assert (tmp_path / "last_summary.json").exists()
    assert (tmp_path / "stride_summary.json").exists()
    assert not Path("md_post_summary.json").exists()
