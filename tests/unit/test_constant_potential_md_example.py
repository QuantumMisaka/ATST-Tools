"""Offline contract checks for the materialized Pt constant-potential MD example."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from ase.constraints import FixAtoms
from ase.io import read
import pytest

from atst_tools.utils.config import ConfigLoader
from atst_tools.calculators.constant_potential import ConstantPotentialCalculator
from atst_tools.calculators.factory import CalculatorFactory


EXAMPLE = Path(__file__).resolve().parents[2] / "examples/19_constant_potential_Pt"


def test_pt_md_materializer_requires_explicit_target(monkeypatch, tmp_path):
    """The example may not silently borrow a target from another Hamiltonian."""
    monkeypatch.syspath_prepend(str(EXAMPLE))
    from materialize_md import materialize_md

    with pytest.raises(TypeError, match="target_mu_ev"):
        materialize_md(tmp_path / "missing")
    with pytest.raises(ValueError, match="finite"):
        materialize_md(tmp_path / "nonfinite", float("nan"))
    assert not (tmp_path / "nonfinite").exists()


def test_pt_md_materializer_preserves_boundary_assets_and_mobile_dofs(monkeypatch, tmp_path):
    """Generated YAML and constrained slab are ready for a site ABACUS run."""
    monkeypatch.syspath_prepend(str(EXAMPLE))
    from materialize_md import materialize_md

    target = 14.125
    destination = tmp_path / "pt_md"
    path = materialize_md(destination, target, steps=3)
    monkeypatch.chdir(destination)
    raw = ConfigLoader.load(str(path))
    normalized = ConfigLoader.normalize(raw)
    calculation = normalized["calculation"]
    calculator = normalized["calculator"]
    parameters = calculator["abacus"]["parameters"]
    cp = calculator["constant_potential"]

    assert calculation["type"] == "md"
    assert calculation["driver"] == "ase"
    assert calculation["ensemble"] == "nvt"
    assert calculation["algorithm"] == "langevin"
    assert calculation["steps"] == 3
    assert cp["energy_boundary"] == "compensated_gate"
    assert cp["reference_electrode"] == "custom"
    assert cp["target_mu_ev"] == target
    assert cp["reference_electrons"] == 216
    assert parameters["ks_solver"] == "cusolver"
    assert parameters["basis_type"] == "lcao"
    assert parameters["gate_flag"] == 1
    assert parameters["dip_cor_flag"] == 1
    assert parameters["block_height"] == 2.0
    assert parameters["out_chg"] == "1 12"
    assert parameters["cal_force"] == 1
    decorated = CalculatorFactory.get_calculator(
        "abacus", normalized, directory=str(destination / "run_md")
    )
    assert isinstance(decorated, ConstantPotentialCalculator)
    assert decorated.target_mu_ev == target

    atoms = read(destination / calculation["init_structure"])
    assert len(atoms) == 12
    fixed = set()
    for constraint in atoms.constraints:
        if isinstance(constraint, FixAtoms):
            fixed.update(constraint.get_indices().tolist())
    assert fixed == set(range(8))
    assert len(atoms) - len(fixed) == 4
    assert atoms.get_volume() > 0

    fixture = destination / "fixture"
    record = json.loads((fixture / "fixture.json").read_text(encoding="utf-8"))
    assert record["target_mu_ev"] == target
    assert record["mobile_atom_indices_1based"] == [9, 10, 11, 12]
    for name, digest in record["sha256"].items():
        assert hashlib.sha256((fixture / name).read_bytes()).hexdigest() == digest
    for mapping in (calculator["abacus"]["pseudopotentials"], calculator["abacus"]["basissets"]):
        for filename in mapping.values():
            assert (fixture / filename).is_file()
