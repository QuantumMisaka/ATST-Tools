"""Ordinary ABACUS MD must not integrate forces without SCF validity."""
import json

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.io import write


@pytest.mark.parametrize("status", [True, False, None])
def test_abacus_md_requires_each_scf_to_converge(tmp_path, monkeypatch, status):
    from atst_tools.workflows import md

    class Electronic(Calculator):
        implemented_properties = ["energy", "forces", "stress"]

        def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            self.results = {"energy": 1., "forces": np.zeros((len(atoms), 3)),
                            "stress": np.zeros(6)}
            if status is not None:
                self.results["scf_converged"] = status

    monkeypatch.chdir(tmp_path)
    atoms = Atoms("H2", positions=[[0, 0, 0], [1, 0, 0]], cell=[8, 8, 8], pbc=True)
    atoms.set_momenta([[.1, 0, 0], [-.1, 0, 0]])
    write("init.traj", atoms)
    monkeypatch.setattr(md.CalculatorFactory, "get_calculator", lambda *a, **kw: Electronic())
    calculation = {"type": "md", "driver": "ase", "ensemble": "nve",
                   "algorithm": "velocityverlet", "init_structure": "init.traj", "steps": 2}
    workflow = md.MDWorkflow({"calculator": {"name": "abacus"}}, "abacus", calculation)
    if status is True:
        workflow.run()
        summary = json.loads((tmp_path / "md_summary.json").read_text())
        assert summary["scf_converged"] is True
    else:
        with pytest.raises(RuntimeError, match="SCF"):
            workflow.run()
        assert not (tmp_path / "md_final.traj").exists()
        assert not (tmp_path / "atst_artifacts.json").exists()
