"""Cross-version contracts for durable constant-potential facts-v1."""

from __future__ import annotations

import copy

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.io import read, write

import cp_v228_reader_fixture as legacy
from atst_tools.calculators import constant_potential as current


class FixedElectronBackend(Calculator):
    """Provide a deterministic converged electronic result for one electron count."""

    implemented_properties = ["energy", "free_energy", "forces"]

    def __init__(self, nelec: float):
        super().__init__()
        self.nelec = nelec

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        energy = 0.5 * self.nelec**2
        self.results.update(
            energy=energy,
            free_energy=energy,
            forces=np.zeros((len(atoms), 3)),
            efermi=3.2 + 0.12 * self.nelec,
            vacuum_level=10.0,
            scf_converged=True,
        )


def _completed_frame(tmp_path):
    """Calculate one converged reference-boundary frame with current code."""
    atoms = Atoms("H", positions=[[0.0, 0.0, 0.0]], cell=[10.0] * 3, pbc=True)
    calculator = current.ConstantPotentialCalculator(
        lambda nelec, directory: FixedElectronBackend(nelec),
        potential_v=1.0,
        reference_electrons=10.0,
        initial_electrons=8.0,
        work_ref=4.6,
        work_ref_source="test",
        capacitance_initial=1.0 / 0.12,
        potential_tolerance_v=1e-12,
        directory=tmp_path,
    )
    atoms.calc = calculator
    atoms.get_potential_energy()
    return atoms, calculator


def test_current_frame_is_readable_by_v228_and_restart_uses_evaluated_count(tmp_path):
    atoms, calculator = _completed_frame(tmp_path)
    trajectory = tmp_path / "checkpoint.traj"
    write(trajectory, atoms)
    restored = read(trajectory)
    expected = calculator.constant_potential_identity

    old_payload = legacy.read_constant_potential_facts(restored)
    assert old_payload is not None
    assert old_payload["facts"]["profile_status"] == "reference_only_unvalidated"
    assert legacy.constant_potential_restart_initial_electrons(
        restored, expected, 1e-12
    ) == pytest.approx(calculator.results["nelec"])


def test_v228_and_transition_frames_remain_readable_by_current_reader(tmp_path):
    atoms, calculator = _completed_frame(tmp_path)
    results = dict(calculator.results)
    results["profile_status"] = "reference_only_unvalidated"
    legacy.publish_constant_potential_facts(
        atoms, results, calculator.constant_potential_identity
    )
    assert current.read_constant_potential_facts(atoms) is not None
    assert current.constant_potential_restart_initial_electrons(
        atoms, calculator.constant_potential_identity, 1e-12
    ) == pytest.approx(results["nelec"])

    malformed = copy.deepcopy(atoms.info[current.CP_FACTS_INFO_KEY])
    malformed["facts"]["profile_status"] = ""
    atoms.info[current.CP_FACTS_INFO_KEY] = malformed
    assert current.read_constant_potential_facts(atoms) is None

    transition = copy.deepcopy(malformed)
    transition["facts"].pop("profile_status")
    atoms.info[current.CP_FACTS_INFO_KEY] = transition
    assert legacy.read_constant_potential_facts(atoms) is None
    assert current.read_constant_potential_facts(atoms) is not None
    assert current.constant_potential_restart_initial_electrons(
        atoms, calculator.constant_potential_identity, 1e-12
    ) == pytest.approx(results["nelec"])


def test_compensated_gate_uses_historical_compatibility_label():
    atoms = Atoms("H", positions=[[0.0, 0.0, 0.0]])
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
    results = {
        "energy": 1.0,
        "free_energy": 1.0,
        "raw_energy": 1.0,
        "raw_free_energy": 1.0,
        "efermi": -15.8,
        "vacuum_level": None,
        "fermishift": None,
        "energy_boundary": "compensated_gate",
        "boundary_parameters": {},
        "compensation": {
            "boundary": "compensated_gate",
            "gate_derivative_ev": 0.4,
            "dipole_derivative_ev": 0.4,
            "compensation_derivative_ev": 0.8,
            "integrated_electrons": 10.0,
            "density_electron_error": 0.0,
            "ionic_valence_electrons": 10.0,
            "total_dipole_ry_au": 0.0,
            "density_precision": 12,
            "density_electron_tolerance": 1e-8,
        },
        "mu_calc": -15.0,
        "mu_target": -15.0,
        "potential_calc": None,
        "potential_target": None,
        "residual_mu": 0.0,
        "residual_v": 0.0,
        "omega_correction": 0.0,
        "nelec": 10.0,
        "reference_electrons": 10.0,
        "delta_nelec": 0.0,
        "forces": np.zeros((1, 3)),
        "scf_converged": True,
        "cp_converged": True,
        "cp_iterations": 1,
        "cp_history": [{"converged": True}],
    }
    payload = current.publish_constant_potential_facts(atoms, results, identity)
    assert payload["facts"]["profile_status"] == (
        "compensated_gate_pending_scientific_acceptance"
    )
    assert legacy.read_constant_potential_facts(atoms) is not None
    assert legacy.constant_potential_restart_initial_electrons(
        atoms, identity, 1e-12
    ) == pytest.approx(10.0)


@pytest.mark.parametrize("reader", [legacy, current], ids=["v228", "current"])
def test_both_readers_reject_stale_geometry_and_identity(tmp_path, reader):
    atoms, calculator = _completed_frame(tmp_path)
    expected = calculator.constant_potential_identity
    assert reader.read_constant_potential_facts(atoms) is not None

    wrong_identity = copy.deepcopy(expected)
    wrong_identity["potential_v"] += 1.0
    with pytest.raises(reader.ConstantPotentialError, match="identity"):
        reader.constant_potential_restart_initial_electrons(atoms, wrong_identity, 1e-12)

    atoms.positions[0, 0] += 0.1
    assert reader.read_constant_potential_facts(atoms) is None
    with pytest.raises(reader.ConstantPotentialError, match="geometry"):
        reader.constant_potential_restart_initial_electrons(atoms, expected, 1e-12)
