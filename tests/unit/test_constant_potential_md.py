"""Runtime and restart contracts for strict constant-potential MD."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.io import read, write

from atst_tools.calculators.constant_potential import (
    ConstantPotentialCalculator,
    ConstantPotentialConvergenceError,
    fixed_hamiltonian_identity_for_config,
    read_constant_potential_facts,
)
from atst_tools.workflows.md import AseMDRunner, AbacusNativeMDRunner


class _ElectronicBackend(Calculator):
    implemented_properties = ["energy", "free_energy", "forces"]

    def __init__(self, nelec: float, *, fail_after_x: float | None = None):
        super().__init__()
        self.nelec = nelec
        self.fail_after_x = fail_after_x

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        x = float(atoms.positions[0, 0])
        q = self.nelec - 10.0 - 0.2 * x
        energy = 0.5 * x**2 + 0.5 * q**2
        forces = np.zeros((len(atoms), 3))
        forces[0, 0] = -x + 0.2 * q
        self.results.update(
            energy=energy,
            free_energy=energy,
            forces=forces,
            efermi=q,
            scf_converged=self.fail_after_x is None or x < self.fail_after_x,
        )


def _setup(monkeypatch, tmp_path, *, steps=3, restart=False, interval=1, fail_after_x=None):
    """Build one MD segment with a real CP loop and analytic inner backend."""
    from atst_tools.workflows import md

    tmp_path.mkdir(parents=True, exist_ok=True)
    monkeypatch.chdir(tmp_path)

    monkeypatch.setattr(
        "atst_tools.utils.gate_compensation.read_gate_compensation",
        lambda directory, nelec, reference_electrons: {
            "boundary": "compensated_gate",
            "gate_derivative_ev": 0.0,
            "dipole_derivative_ev": 0.0,
            "compensation_derivative_ev": 0.0,
            "integrated_electrons": nelec,
            "density_electron_error": 0.0,
            "ionic_valence_electrons": reference_electrons,
            "total_dipole_ry_au": 0.0,
            "density_precision": 12,
            "density_electron_tolerance": 1e-8,
        },
    )
    run_dir = tmp_path / "run"
    initial = tmp_path / "initial.traj"
    if not initial.exists():
        atoms = Atoms("H", positions=[[0.5, 0, 0]], cell=[10, 10, 10], pbc=True)
        atoms.set_momenta([[0.1, 0, 0]])
        write(initial, atoms)
    cp = {
        "energy_boundary": "compensated_gate",
        "reference_electrode": "custom",
        "target_mu_ev": 1.0,
        "reference_electrons": 10.0,
        "potential_tolerance_v": 1e-9,
    }
    config = {"calculator": {"name": "abacus", "constant_potential": cp}}
    calc_config = {
        "type": "md",
        "driver": "ase",
        "ensemble": "nve",
        "algorithm": "velocityverlet",
        "init_structure": str(initial),
        "steps": steps,
        "timestep_fs": 0.1,
        "temperature_K": 300.0,
        "directory": str(run_dir),
        "trajectory": str(run_dir / "md.traj"),
        "logfile": str(run_dir / "md.log"),
        "loginterval": interval,
        "summary_file": str(run_dir / "md_summary.json"),
        "final_structure": str(run_dir / "md_final.traj"),
        "artifact_manifest": str(run_dir / "atst_artifacts.json"),
        "postprocess": {"summary": {"enabled": False}, "convert": {"enabled": False}},
        "restart": restart,
    }
    calculators = []

    def make_calculator(*args, **kwargs):
        calculator = ConstantPotentialCalculator(
            lambda nelec, directory: _ElectronicBackend(nelec, fail_after_x=fail_after_x),
            potential_v=None,
            reference_electrons=10.0,
            initial_electrons=10.0,
            energy_boundary="compensated_gate",
            target_mu_ev=1.0,
            reference_electrode="custom",
            capacitance_initial=1.0,
            capacitance_unit="e/V",
            potential_tolerance_v=1e-9,
            fixed_hamiltonian_identity=fixed_hamiltonian_identity_for_config(config),
            directory=run_dir,
        )
        calculators.append(calculator)
        return calculator

    monkeypatch.setattr(md.CalculatorFactory, "get_calculator", make_calculator)
    return AseMDRunner(config, "abacus", calc_config), calculators


def test_cp_md_writes_validated_steps_and_restart_keeps_global_cadence(monkeypatch, tmp_path):
    runner, calculators = _setup(monkeypatch, tmp_path, steps=7, interval=5)
    runner.run()
    root = tmp_path / "run"
    checkpoint = read(root / "cp_md_checkpoint.traj")
    assert checkpoint.info["atst_md_step"] == 7
    assert read_constant_potential_facts(checkpoint) is not None
    assert [frame.info["atst_md_step"] for frame in read(root / "md.traj", index=":")] == [0, 5]
    assert len((root / "cp_md_diagnostics.jsonl").read_text().splitlines()) == 8
    assert calculators[0].initial_electrons != 10.0

    resumed, _ = _setup(monkeypatch, tmp_path, steps=4, restart=True, interval=5)
    resumed.run()
    frames = read(root / "md.traj", index=":")
    assert [frame.info["atst_md_step"] for frame in frames] == [0, 5, 10]
    assert read(root / "cp_md_checkpoint.traj").info["atst_md_step"] == 11
    assert len((root / "cp_md_diagnostics.jsonl").read_text().splitlines()) == 12
    summary = json.loads((root / "md_summary.json").read_text())
    assert summary["completed_steps"] == 11


def test_cp_md_rejects_bad_restart_before_modifying_outputs(monkeypatch, tmp_path):
    runner, _ = _setup(monkeypatch, tmp_path, steps=1)
    runner.run()
    root = tmp_path / "run"
    checkpoint_file = root / "cp_md_checkpoint.traj"
    checkpoint = read(checkpoint_file)
    checkpoint.arrays.pop("momenta")
    write(checkpoint_file, checkpoint)
    old_trajectory = (root / "md.traj").read_bytes()
    resumed, _ = _setup(monkeypatch, tmp_path, steps=1, restart=True)
    with pytest.raises(ValueError, match="momenta"):
        resumed.run()
    assert (root / "md.traj").read_bytes() == old_trajectory


def test_cp_md_first_step_failure_retains_only_step_zero(monkeypatch, tmp_path):
    runner, _ = _setup(monkeypatch, tmp_path, steps=2, fail_after_x=0.50001)
    with pytest.raises((ConstantPotentialConvergenceError, RuntimeError)):
        runner.run()
    root = tmp_path / "run"
    assert read(root / "cp_md_checkpoint.traj").info["atst_md_step"] == 0
    assert [frame.info["atst_md_step"] for frame in read(root / "md.traj", index=":")] == [0]
    assert json.loads((root / "md_summary.json").read_text())["status"] == "failed"


def test_cp_md_fresh_run_rejects_prior_outputs_without_modifying_them(monkeypatch, tmp_path):
    runner, _ = _setup(monkeypatch, tmp_path, steps=1)
    runner.run()
    root = tmp_path / "run"
    before = {name: (root / name).read_bytes() for name in (
        "md.traj", "md.log", "cp_md_checkpoint.traj", "cp_md_diagnostics.jsonl"
    )}
    again, _ = _setup(monkeypatch, tmp_path, steps=1)
    with pytest.raises(ValueError, match="overwrite existing outputs"):
        again.run()
    assert {name: (root / name).read_bytes() for name in before} == before


def test_cp_md_rejects_same_step_foreign_trajectory_tail(monkeypatch, tmp_path):
    runner, _ = _setup(monkeypatch, tmp_path / "first", steps=1)
    runner.run()
    other, _ = _setup(monkeypatch, tmp_path / "second", steps=1)
    initial = read(tmp_path / "second" / "initial.traj")
    initial.positions[0, 0] += 0.2
    write(tmp_path / "second" / "initial.traj", initial)
    other.run()
    first_root = tmp_path / "first" / "run"
    # A valid same-step frame from another geometry must not be accepted.
    original = read(first_root / "md.traj", index=":")
    replacement = read(tmp_path / "second" / "run" / "md.traj", index=-1)
    write(first_root / "md.traj", [original[0], replacement])
    before = (first_root / "md.traj").read_bytes()
    resumed, _ = _setup(monkeypatch, tmp_path / "first", steps=1, restart=True)
    with pytest.raises(ValueError, match="geometry disagrees"):
        resumed.run()
    assert (first_root / "md.traj").read_bytes() == before


def test_cp_md_rejects_log_tail_conflict_before_append(monkeypatch, tmp_path):
    runner, _ = _setup(monkeypatch, tmp_path, steps=2)
    runner.run()
    root = tmp_path / "run"
    lines = (root / "md.log").read_text().splitlines()
    fields = lines[-1].split()
    fields[6] = str(float(fields[6]) + 1.0)
    lines[-1] = " ".join(fields)
    (root / "md.log").write_text("\n".join(lines) + "\n")
    before = (root / "md.log").read_bytes()
    resumed, _ = _setup(monkeypatch, tmp_path, steps=1, restart=True)
    with pytest.raises(ValueError, match="log tail disagrees"):
        resumed.run()
    assert (root / "md.log").read_bytes() == before


def test_cp_md_rejects_foreign_regular_tail_before_unscheduled_checkpoint(monkeypatch, tmp_path):
    first, _ = _setup(monkeypatch, tmp_path / "first", steps=7, interval=5)
    first.run()
    second, _ = _setup(monkeypatch, tmp_path / "second", steps=7, interval=5)
    initial = read(tmp_path / "second" / "initial.traj")
    initial.positions[0, 1] += 0.2  # The analytic electronic model depends only on x.
    write(tmp_path / "second" / "initial.traj", initial)
    second.run()
    root = tmp_path / "first" / "run"
    own_frames = read(root / "md.traj", index=":")
    foreign_tail = read(tmp_path / "second" / "run" / "md.traj", index=-1)
    assert own_frames[-1].info["atst_md_step"] == 5
    assert foreign_tail.info["atst_md_step"] == 5
    write(root / "md.traj", [own_frames[0], foreign_tail])
    previous = (root / "md.traj").read_bytes()

    resumed, _ = _setup(monkeypatch, tmp_path / "first", steps=1, restart=True, interval=5)
    with pytest.raises(ValueError, match="nuclear state"):
        resumed.run()
    assert (root / "md.traj").read_bytes() == previous


def test_cp_md_rejects_tampered_checkpoint_grand_potential_row(monkeypatch, tmp_path):
    runner, _ = _setup(monkeypatch, tmp_path, steps=7, interval=5)
    runner.run()
    root = tmp_path / "run"
    diagnostic = root / "cp_md_diagnostics.jsonl"
    rows = [json.loads(line) for line in diagnostic.read_text().splitlines()]
    rows[-1]["grand_potential_eV"] = 123456.0
    diagnostic.write_text("\n".join(json.dumps(row) for row in rows) + "\n")
    before = diagnostic.read_bytes()

    resumed, _ = _setup(monkeypatch, tmp_path, steps=1, restart=True, interval=5)
    with pytest.raises(ValueError, match="grand potential"):
        resumed.run()
    assert diagnostic.read_bytes() == before


@pytest.mark.parametrize("driver", [None, "ase", "abacus_native"])
def test_native_runner_rejects_cp_before_any_input_preparation(monkeypatch, tmp_path, driver):
    """A direct native runner must not trust the caller's driver field."""
    calculation = {"directory": str(tmp_path / "native")}
    if driver is not None:
        calculation["driver"] = driver
    runner = AbacusNativeMDRunner(
        {
            "calculator": {
                "name": "abacus",
                "constant_potential": {
                    "energy_boundary": "compensated_gate",
                    "reference_electrode": "custom",
                    "target_mu_ev": 1.0,
                },
            }
        },
        "abacus",
        calculation,
    )
    monkeypatch.setattr(runner, "_prepare_inputs", lambda: pytest.fail("native inputs prepared"))
    monkeypatch.setattr(runner, "_start_process", lambda: pytest.fail("native engine started"))

    with pytest.raises(ValueError, match="constant-potential MD.*driver=ase"):
        runner.run()
    assert not (tmp_path / "native").exists()
