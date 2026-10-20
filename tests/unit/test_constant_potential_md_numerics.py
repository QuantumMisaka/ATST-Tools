"""Independent analytic checks of CP forces during real ASE MD integration."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.constraints import FixAtoms, Hookean
from ase.io import read, write

from atst_tools.calculators.constant_potential import (
    ConstantPotentialCalculator,
    ConstantPotentialConvergenceError,
    ConstantPotentialEvaluationError,
    fixed_hamiltonian_identity_for_config,
    read_constant_potential_facts,
)
from atst_tools.workflows import md


K = 2.0
C = 0.4
A = 0.3
N0 = 10.0
MU = 0.2
GATE_DERIVATIVE = 0.35


def stationary_electrons(atoms: Atoms) -> float:
    """Return the independent electronic minimizer at one geometry."""
    return N0 + A * atoms.positions[1, 0] + C * MU


def minimized_omega(atoms: Atoms) -> float:
    """Return the closed-form grand potential of the analytic Hamiltonian."""
    return 0.5 * K * float(np.sum(atoms.positions**2)) - MU * A * atoms.positions[1, 0] - 0.5 * C * MU**2


def stationary_forces(atoms: Atoms) -> np.ndarray:
    """Return the exact force after minimizing the electron number."""
    forces = -K * atoms.positions.copy()
    forces[1, 0] += A * MU
    return forces


class AnalyticElectronicBackend(Calculator):
    """Fixed-electron backend with coupled geometry and electron response."""

    implemented_properties = ["energy", "free_energy", "forces"]

    def __init__(self, nelec: float, directory: str, trace: list[tuple[tuple[float, ...], float]], *, fail: bool = False):
        super().__init__(directory=directory)
        self.nelec = float(nelec)
        self.trace = trace
        self.fail = fail

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        geometry = tuple(np.asarray(atoms.positions).ravel())
        self.trace.append((geometry, self.nelec))
        displacement = self.nelec - N0 - A * atoms.positions[1, 0]
        energy = 0.5 * K * float(np.sum(atoms.positions**2)) + displacement**2 / (2 * C)
        forces = -K * atoms.positions.copy()
        forces[1, 0] += A * displacement / C
        self.results.update(
            energy=energy,
            free_energy=energy,
            forces=forces,
            efermi=displacement / C - GATE_DERIVATIVE,
            scf_converged=not self.fail,
        )


def compensation(directory: str, nelec: float, reference_electrons: float) -> dict:
    """Supply an internally consistent synthetic gate-density measurement."""
    return {
        "boundary": "compensated_gate",
        "gate_derivative_ev": GATE_DERIVATIVE,
        "dipole_derivative_ev": 0.0,
        "compensation_derivative_ev": GATE_DERIVATIVE,
        "integrated_electrons": float(nelec),
        "density_electron_error": 0.0,
        "ionic_valence_electrons": float(reference_electrons),
        "total_dipole_ry_au": 0.0,
        "density_precision": 12,
        "density_electron_tolerance": 1e-8,
    }


def cp_calculator(
    root: Path,
    trace: list,
    *,
    fail: bool = False,
    max_iterations: int = 8,
    config: dict | None = None,
) -> ConstantPotentialCalculator:
    """Build the production CP solver around the analytic fixed-electron backend."""
    return ConstantPotentialCalculator(
        lambda nelec, directory: AnalyticElectronicBackend(nelec, directory, trace, fail=fail),
        potential_v=None,
        energy_boundary="compensated_gate",
        target_mu_ev=MU,
        reference_electrode="custom",
        reference_electrons=N0,
        initial_electrons=N0,
        potential_tolerance_v=1e-10,
        capacitance_initial=C,
        capacitance_unit="e/V",
        max_iterations=max_iterations,
        directory=root,
        fixed_hamiltonian_identity=(
            fixed_hamiltonian_identity_for_config(config) if config is not None else None
        ),
    )


def initial_atoms(*, constrained: bool = True) -> Atoms:
    """Create a small fixed-cell system with mobile nuclear degrees of freedom."""
    atoms = Atoms(
        "H3",
        positions=[[0.0, 0.0, 0.0], [0.2, 0.08, -0.04], [-0.12, 0.05, 0.02]],
        cell=[10.0, 10.0, 10.0],
        pbc=True,
    )
    if constrained:
        atoms.set_constraint(FixAtoms(indices=[0]))
    atoms.set_momenta([[0.0, 0.0, 0.0], [0.12, -0.07, 0.03], [-0.08, 0.05, -0.02]])
    return atoms


def md_config(tmp_path: Path, *, algorithm: str, steps: int = 4, timestep_fs: float = 0.25) -> tuple[dict, dict]:
    """Build the direct workflow config while retaining the real CP decorator."""
    (tmp_path / "run").mkdir(exist_ok=True)
    calculation = {
        "type": "md",
        "driver": "ase",
        "init_structure": str(tmp_path / "init.traj"),
        "directory": str(tmp_path / "run"),
        "ensemble": "nve" if algorithm == "velocityverlet" else "nvt",
        "algorithm": algorithm,
        "steps": steps,
        "timestep_fs": timestep_fs,
        "temperature_K": 100.0,
        "stationary": False,
        "friction_fs_inv": 0.01,
        "taut_fs": 10.0,
        "seed": 5,
        "trajectory": str(tmp_path / "run" / "md.traj"),
        "logfile": str(tmp_path / "run" / "md.log"),
        "summary_file": str(tmp_path / "run" / "md_summary.json"),
        "final_structure": str(tmp_path / "run" / "md_final.traj"),
        "artifact_manifest": str(tmp_path / "run" / "atst_artifacts.json"),
    }
    config = {
        "calculation": calculation,
        "calculator": {
            "name": "abacus",
            "abacus": {"parameters": {"nelec": N0}},
            "constant_potential": {
                "energy_boundary": "compensated_gate",
                "reference_electrode": "custom",
                "target_mu_ev": MU,
                "reference_electrons": N0,
                "potential_tolerance_v": 1e-10,
            },
        },
    }
    return config, calculation


def test_compensated_cp_force_matches_minimized_grand_potential_difference(monkeypatch, tmp_path):
    """The force must be a derivative of the minimized energy at each geometry."""
    monkeypatch.setattr("atst_tools.utils.gate_compensation.read_gate_compensation", compensation)
    samples = []
    energies = []
    for shift in (-1e-3, 0.0, 1e-3):
        atoms = initial_atoms()
        atoms.positions[1, 0] += shift
        atoms.calc = cp_calculator(tmp_path / f"x{shift:+.3f}", samples)
        omega = atoms.get_potential_energy(force_consistent=True)
        energies.append(omega)
        assert omega == pytest.approx(minimized_omega(atoms), abs=2e-12)
        assert atoms.calc.results["nelec"] == pytest.approx(stationary_electrons(atoms), abs=2e-11)
        if shift == 0.0:
            force = atoms.get_forces()[1, 0]
            assert force == pytest.approx(stationary_forces(atoms)[1, 0], abs=2e-11)
    # The centered difference uses two independently reconverged CP solves.
    center_difference = -(energies[2] - energies[0]) / 2e-3
    assert force == pytest.approx(center_difference, abs=1e-10)


def test_nonconverged_electronic_state_publishes_no_force(monkeypatch, tmp_path):
    """One insufficient CP iteration must fail before MD can consume a force."""
    monkeypatch.setattr("atst_tools.utils.gate_compensation.read_gate_compensation", compensation)
    atoms = initial_atoms()
    atoms.calc = cp_calculator(tmp_path / "nonconverged", [], max_iterations=1)
    with pytest.raises(ConstantPotentialConvergenceError):
        atoms.get_forces()
    assert "forces" not in atoms.calc.results
    assert read_constant_potential_facts(atoms) is None


def test_failed_scf_state_publishes_no_force(monkeypatch, tmp_path):
    """A reported SCF failure cannot become a force even with a valid mu."""
    monkeypatch.setattr("atst_tools.utils.gate_compensation.read_gate_compensation", compensation)
    atoms = initial_atoms()
    atoms.calc = cp_calculator(tmp_path / "failed_scf", [], fail=True)
    with pytest.raises(ConstantPotentialEvaluationError, match="not converged"):
        atoms.get_forces()
    assert "forces" not in atoms.calc.results
    assert read_constant_potential_facts(atoms) is None


def test_hookean_constraint_contributes_once_to_energy_and_force(monkeypatch, tmp_path):
    """The private electronic probe must not duplicate nuclear constraint work."""
    monkeypatch.setattr("atst_tools.utils.gate_compensation.read_gate_compensation", compensation)
    atoms = initial_atoms(constrained=False)
    atoms.set_constraint(Hookean(a1=1, a2=2, rt=0.1, k=3.0))
    atoms.calc = cp_calculator(tmp_path / "hookean", [])
    displacement = atoms.positions[2] - atoms.positions[1]
    distance = float(np.linalg.norm(displacement))
    extension = distance - 0.1
    hook_energy = 0.5 * 3.0 * extension**2
    expected_forces = stationary_forces(atoms)
    hook_force = 3.0 * extension * displacement / distance
    expected_forces[1] += hook_force
    expected_forces[2] -= hook_force
    assert atoms.get_potential_energy(force_consistent=True) == pytest.approx(
        minimized_omega(atoms) + hook_energy, abs=2e-11
    )
    np.testing.assert_allclose(atoms.get_forces(), expected_forces, atol=2e-11)


@pytest.mark.parametrize("algorithm", ["velocityverlet", "bussi", "langevin", "nvtberendsen"])
def test_md_frames_follow_analytic_cp_surface(monkeypatch, tmp_path, algorithm):
    """Real ASE integrators consume only same-geometry converged CP forces."""
    monkeypatch.setattr("atst_tools.utils.gate_compensation.read_gate_compensation", compensation)
    trace: list[tuple[tuple[float, ...], float]] = []
    calculators = []

    def calculator_factory(*args, **kwargs):
        result = cp_calculator(tmp_path / "run", trace)
        calculators.append(result)
        return result

    monkeypatch.setattr(md.CalculatorFactory, "get_calculator", calculator_factory)
    write(tmp_path / "init.traj", initial_atoms())
    config, calculation = md_config(tmp_path, algorithm=algorithm)
    monkeypatch.chdir(tmp_path)
    md.MDWorkflow(config, "abacus", calculation).run()

    frames = read(calculation["trajectory"], index=":")
    rows = [json.loads(line) for line in (tmp_path / "run" / "cp_md_diagnostics.jsonl").read_text().splitlines()]
    assert [row["step"] for row in rows] == list(range(5))
    assert len(frames) == 5
    assert len(calculators) == 1
    groups: list[tuple[tuple[float, ...], list[float]]] = []
    for geometry, nelec in trace:
        if not groups or groups[-1][0] != geometry:
            groups.append((geometry, []))
        groups[-1][1].append(nelec)
    # The backend calls at each geometry must match the solver's declared
    # iteration count. Logger/trajectory reads must not start another batch.
    assert len(groups) == 5
    assert [len(samples) for _, samples in groups] == [row["cp_iterations"] for row in rows]
    assert all(row["cp_iterations"] >= 1 for row in rows)
    for previous, current in zip(rows, groups[1:]):
        assert current[1][0] == pytest.approx(previous["nelec"], abs=2e-11)
    for step, (frame, row) in enumerate(zip(frames, rows)):
        assert frame.info["atst_md_step"] == step
        assert frame.info["atst_md_time_fs"] == pytest.approx(step * calculation["timestep_fs"])
        facts = read_constant_potential_facts(frame)
        assert facts is not None
        measured = facts["facts"]
        assert measured["nelec"] == pytest.approx(stationary_electrons(frame), abs=2e-9)
        assert measured["residual_mu"] == pytest.approx(0.0, abs=1e-10)
        np.testing.assert_allclose(measured["forces"], stationary_forces(frame), atol=2e-9)
        assert measured["free_energy"] == pytest.approx(minimized_omega(frame), abs=2e-9)
        assert row["nelec"] == pytest.approx(measured["nelec"])
        assert row["grand_potential_eV"] == pytest.approx(measured["free_energy"])
        np.testing.assert_allclose(frame.positions[0], initial_atoms().positions[0], atol=0.0)
    checkpoint = read(tmp_path / "run" / "cp_md_checkpoint.traj")
    assert checkpoint.info["atst_md_step"] == 4
    assert read_constant_potential_facts(checkpoint) is not None
    summary = json.loads(Path(calculation["summary_file"]).read_text())
    assert summary["completed_steps"] == 4
    manifest = json.loads(Path(calculation["artifact_manifest"]).read_text())
    assert {item["role"] for item in manifest["artifacts"]} >= {"cp_md_checkpoint", "cp_md_diagnostics"}


def test_verlet_energy_error_shrinks_with_timestep(monkeypatch, tmp_path):
    """Second-order Verlet refines conservation of K plus grand potential."""
    monkeypatch.setattr("atst_tools.utils.gate_compensation.read_gate_compensation", compensation)
    errors = []
    for label, timestep, steps in (("coarse", 0.5, 40), ("fine", 0.25, 80)):
        root = tmp_path / label
        root.mkdir()
        write(root / "init.traj", initial_atoms(constrained=False))
        monkeypatch.setattr(
            md.CalculatorFactory,
            "get_calculator",
            lambda *args, _root=root, **kwargs: cp_calculator(_root / "run", []),
        )
        config, calculation = md_config(root, algorithm="velocityverlet", steps=steps, timestep_fs=timestep)
        monkeypatch.chdir(root)
        md.MDWorkflow(config, "abacus", calculation).run()
        rows = [json.loads(line) for line in (root / "run" / "cp_md_diagnostics.jsonl").read_text().splitlines()]
        energies = np.array([row["grand_potential_eV"] + row["kinetic_energy_eV"] for row in rows])
        errors.append(float(np.max(np.abs(energies - energies[0]))))
    assert errors[0] > 1e-8
    assert errors[1] < 0.4 * errors[0]


def test_bussi_with_one_mobile_atom_and_three_free_dofs(monkeypatch, tmp_path):
    """A constrained slab still has a valid thermostat when only one atom moves."""
    monkeypatch.setattr("atst_tools.utils.gate_compensation.read_gate_compensation", compensation)
    atoms = initial_atoms(constrained=False)
    atoms.set_constraint(FixAtoms(indices=[0, 2]))
    atoms.set_momenta([[0.0, 0.0, 0.0], [0.12, -0.07, 0.03], [0.0, 0.0, 0.0]])
    write(tmp_path / "init.traj", atoms)
    monkeypatch.setattr(
        md.CalculatorFactory,
        "get_calculator",
        lambda *args, **kwargs: cp_calculator(tmp_path / "run", []),
    )
    config, calculation = md_config(tmp_path, algorithm="bussi", steps=2)
    monkeypatch.chdir(tmp_path)
    md.MDWorkflow(config, "abacus", calculation).run()
    frames = read(calculation["trajectory"], index=":")
    assert [frame.info["atst_md_step"] for frame in frames] == [0, 1, 2]
    for frame in frames:
        np.testing.assert_allclose(frame.positions[[0, 2]], atoms.positions[[0, 2]], atol=0.0)
        assert read_constant_potential_facts(frame) is not None


def test_analytic_verlet_restart_matches_continuous_global_cadence(monkeypatch, tmp_path):
    """A split at an unscheduled step preserves state and regular output phase."""
    monkeypatch.setattr("atst_tools.utils.gate_compensation.read_gate_compensation", compensation)
    results = {}
    for name in ("continuous", "split"):
        root = tmp_path / name
        root.mkdir()
        write(root / "init.traj", initial_atoms())
        monkeypatch.setattr(
            md.CalculatorFactory,
            "get_calculator",
            lambda *args, _root=root, **kwargs: cp_calculator(_root / "run", [], config=config),
        )
        config, calculation = md_config(
            root, algorithm="velocityverlet", steps=8 if name == "continuous" else 3
        )
        calculation["loginterval"] = 2
        monkeypatch.chdir(root)
        md.MDWorkflow(config, "abacus", calculation).run()
        if name == "split":
            checkpoint = read(root / "run" / "cp_md_checkpoint.traj")
            assert checkpoint.info["atst_md_step"] == 3
            assert read(root / "run" / "md.traj", index=-1).info["atst_md_step"] == 2
            calculation["restart"] = True
            calculation["steps"] = 5
            md.MDWorkflow(config, "abacus", calculation).run()
        frames = read(root / "run" / "md.traj", index=":")
        rows = [
            json.loads(line)
            for line in (root / "run" / "cp_md_diagnostics.jsonl").read_text().splitlines()
        ]
        assert [frame.info["atst_md_step"] for frame in frames] == [0, 2, 4, 6, 8]
        assert [row["step"] for row in rows] == list(range(9))
        results[name] = (read(root / "run" / "md_final.traj"), rows)
    full, full_rows = results["continuous"]
    resumed, resumed_rows = results["split"]
    np.testing.assert_allclose(resumed.positions, full.positions, atol=2e-11)
    np.testing.assert_allclose(resumed.get_momenta(), full.get_momenta(), atol=2e-11)
    assert resumed_rows[-1]["nelec"] == pytest.approx(full_rows[-1]["nelec"], abs=2e-11)
    assert resumed_rows[-1]["grand_potential_eV"] == pytest.approx(
        full_rows[-1]["grand_potential_eV"], abs=2e-11
    )
