"""Compare the decorator with a trace executed by the pinned upstream FCP v2."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes

from atst_tools.calculators.constant_potential import ConstantPotentialCalculator

REFERENCE = Path(__file__).parents[1] / "data/constant_potential_reference.json"


class AnalyticAbacus(Calculator):
    """Expose the same analytic response through an ABACUS-style log."""

    implemented_properties = ["energy", "free_energy", "forces"]

    def __init__(self, nelec, directory, trace):
        super().__init__(directory=directory)
        self.nelec = nelec
        self.trace = trace

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        x = self.nelec - 10
        energy = 100 + 5*x + .2*x*x + .005*x**4
        ef = 5 + .4*x + .02*x**3
        self.results.update(energy=energy, free_energy=energy,
                            forces=np.zeros((len(atoms), 3)), efermi=ef,
                            scf_converged=True)
        output = Path(self.directory) / "OUT.ABACUS"
        output.mkdir(parents=True, exist_ok=True)
        (output / "running_scf.log").write_text(
            " charge density convergence is achieved\n"
            f" EFERMI = {ef:.14g} eV\nThe vacuum level is 10.0 eV\n"
        )
        self.trace.append(self.nelec)


def test_fcp_v2_reference_trajectory_and_evaluated_state(tmp_path):
    """Compare the normal FCP v2 trajectory field by field.

    The committed reference uses a deterministic analytic backend, so this is
    an algorithm-level comparison, not a real ase-abacus/ABACUS engine run.
    FCP v2's ``always_adjust=True`` run is retained for provenance only; ATST
    currently stops at the converged state and does not expose that flag.
    """
    reference = json.loads(REFERENCE.read_text())
    trace = []

    def factory(nelec, directory):
        return AnalyticAbacus(nelec, directory, trace)

    calculator = ConstantPotentialCalculator(
        factory, directory=tmp_path, work_ref_source="pinned FCP v2 comparison",
        **reference["parameters"],
    )
    atoms = Atoms("H", positions=[[0, 0, 0]], cell=reference["cell"], pbc=True)
    atoms.calc = calculator
    energy = atoms.get_potential_energy()
    normal, adjusted = reference["runs"]
    expected = normal["trace"]
    target_mu = -(
        reference["parameters"]["potential_v"] + reference["parameters"]["work_ref"]
    )
    # Upstream fits six-decimal text samples; the port retains full precision.
    # The tolerances cover that declared representation change, not a new solver.
    assert len(trace) == len(expected)
    assert len(calculator.evaluation_history) == len(expected)
    for actual, expected_row in zip(calculator.evaluation_history, expected):
        assert actual["nelec"] == pytest.approx(expected_row["nelec"], abs=2e-5)
        assert actual["efermi"] == pytest.approx(expected_row["efermi"], abs=1e-5)
        assert actual["mu_calc"] == pytest.approx(expected_row["mu"], abs=1e-5)
        assert actual["raw_energy"] == pytest.approx(
            expected_row["raw_energy"], abs=1e-4
        )
        assert actual["residual_mu"] == pytest.approx(
            expected_row["mu"] - target_mu, abs=1e-5
        )
    assert calculator.results["nelec"] == pytest.approx(trace[-1], abs=1e-12)
    assert energy == pytest.approx(normal["published_energy"], abs=1e-7)
    assert calculator.results["free_energy"] == pytest.approx(
        normal["published_energy"], abs=1e-7
    )
    assert calculator.results["raw_energy"] == pytest.approx(
        expected[-1]["raw_energy"], abs=1e-4
    )
    assert calculator.results["mu_target"] == pytest.approx(target_mu, abs=1e-12)
    assert calculator.results["mu_calc"] == pytest.approx(expected[-1]["mu"], abs=1e-5)
    assert calculator.results["residual_mu"] == pytest.approx(
        expected[-1]["mu"] - target_mu, abs=1e-5
    )
    assert calculator.results["vacuum_level"] == pytest.approx(10.0, abs=1e-12)
    assert calculator.results["fermishift"] == pytest.approx(-10.0, abs=1e-12)
    assert calculator.results["cp_iterations"] == len(expected)
    assert calculator.results["cp_converged"] is True
    assert adjusted["last_evaluated_nelec"] == normal["last_evaluated_nelec"]
    assert adjusted["next_guess_nelec"] != adjusted["last_evaluated_nelec"]
    count = len(trace)
    np.testing.assert_allclose(atoms.get_forces(), 0)
    assert len(trace) == count


def test_fcp_v2_reference_provenance_is_pinned():
    """Keep the reference tied to the reviewed upstream FCP v2 source."""
    reference = json.loads(REFERENCE.read_text())

    assert reference["reference_revision"] == "9f96cacb6b69cb39aeb180e73f415f806e16a75e"
    assert reference["reference_path"] == "version2/FCPelectrochem.py"
    assert reference["source_sha256"] == (
        "bd1cbf41abf34d68a41f6f1579275e6c082e8f0f8be09c2d629153c2aa745040"
    )
    assert [run["always_adjust"] for run in reference["runs"]] == [False, True]
