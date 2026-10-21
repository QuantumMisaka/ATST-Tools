"""Shared endpoint relaxation stage for NEB-family workflows."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import Any

from atst_tools.utils.convergence import StageRecord, as_finite_float, as_step_count
from atst_tools.utils.neb_endpoints import (
    ENDPOINT_OPTIMIZED,
    freeze_current_results,
    has_compatible_cp_endpoint_results,
)


def relax_neb_endpoints(
    images: Sequence,
    endpoint_config: dict[str, Any],
    get_calculator: Callable[[str, Any], Any],
    optimizer_class: Any,
    *,
    constant_potential_identity: dict[str, Any] | None = None,
) -> list[StageRecord]:
    """Relax missing-result endpoints and preserve optimizer-owned facts.

    The same energy/force and CP energy-definition test determines whether an
    endpoint can skip Relax and whether NEB can reuse it afterward. Existing
    result provenance and calculator parameter identity are not admission rules.
    """
    if not endpoint_config.get("enabled"):
        return []
    skip_existing = endpoint_config.get("skip_if_has_results", True)
    endpoint_fmax = endpoint_config.get("fmax", 0.05)
    endpoint_steps = endpoint_config.get("max_steps", 100)
    stages: list[StageRecord] = []
    for label, atoms in (("initial", images[0]), ("final", images[-1])):
        if skip_existing and has_compatible_cp_endpoint_results(atoms, constant_potential_identity):
            stages.append(
                StageRecord(
                    name=f"endpoint_{label}_relax",
                    role="endpoint",
                    criterion="ase_optimizer",
                    status="skipped",
                )
            )
            continue
        atoms.calc = get_calculator(label, atoms)
        opt = optimizer_class(atoms, trajectory=f"endpoint_{label}_relax.traj")
        stage_converged = opt.run(fmax=endpoint_fmax, steps=endpoint_steps)
        freeze_current_results(atoms, status=ENDPOINT_OPTIMIZED)
        stages.append(
            StageRecord(
                name=f"endpoint_{label}_relax",
                role="endpoint",
                criterion="ase_optimizer",
                converged=stage_converged,
                fmax=as_finite_float(endpoint_fmax),
                steps=as_step_count(endpoint_steps),
                actual_steps=as_step_count(getattr(opt, "nsteps", None)),
            )
        )
    return stages
