"""Regression tests for DMF ordered endpoint identity."""

import pytest
from ase import Atoms
from ase.io import write

from atst_tools.workflows.dmf import DMFWorkflow


@pytest.mark.parametrize(
    ("initial", "final", "message"),
    [
        ("HO", "OH", "ordered chemical symbols"),
        ("HO", "HOH", "same number of atoms"),
    ],
)
def test_dmf_rejects_mismatched_endpoints_before_optimizer_load(
    tmp_path, monkeypatch, initial, final, message
):
    """Reject correspondence errors before constructing an optimization path."""
    init_path = tmp_path / "initial.traj"
    final_path = tmp_path / "final.traj"
    write(init_path, Atoms(initial))
    write(final_path, Atoms(final))
    workflow = DMFWorkflow(
        {}, "dp", {"type": "dmf", "init_file": str(init_path), "final_file": str(final_path)}
    )
    optimizer_loads = []
    monkeypatch.setattr(workflow, "_load_pydmf", lambda: optimizer_loads.append(True))

    with pytest.raises(ValueError, match=message):
        workflow.run()
    assert optimizer_loads == []


def test_dmf_matching_endpoints_reach_optimizer_load(tmp_path, monkeypatch):
    """A valid ordered endpoint pair passes the correspondence guard."""
    init_path = tmp_path / "initial.traj"
    final_path = tmp_path / "final.traj"
    write(init_path, Atoms("HO", positions=[[0, 0, 0], [1, 0, 0]]))
    write(final_path, Atoms("HO", positions=[[0, 0, 0], [2, 0, 0]]))
    workflow = DMFWorkflow(
        {}, "dp", {"type": "dmf", "init_file": str(init_path), "final_file": str(final_path)}
    )

    class OptimizerReached(Exception):
        """Mark the valid boundary without solving IPOPT."""

    def mark_optimizer_load():
        raise OptimizerReached

    monkeypatch.setattr(workflow, "_load_pydmf", mark_optimizer_load)
    with pytest.raises(OptimizerReached):
        workflow.run()
