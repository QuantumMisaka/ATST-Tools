"""Endpoint result validation and preparation for NEB-family workflows."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import Any

import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator

from atst_tools.calculators.constant_potential import (
    CP_FACTS_INFO_KEY,
    CP_IDENTITY_INFO_KEY,
    _mapping_contains,
    clear_constant_potential_facts,
    constant_potential_identity_for_config,
    publish_constant_potential_facts,
    read_constant_potential_facts,
)


ENDPOINT_RESULT_KEY = "atst_endpoint_result"
ENDPOINT_PLACEHOLDER = "placeholder"
ENDPOINT_PROVIDED = "provided"
ENDPOINT_COMPUTED = "computed"
ENDPOINT_OPTIMIZED = "optimized"
CP_ENDPOINT_IDENTITY_KEY = CP_IDENTITY_INFO_KEY


def mark_endpoint_result(atoms, status: str) -> None:
    """Mark how endpoint energy/force results were obtained."""
    atoms.info[ENDPOINT_RESULT_KEY] = status


def is_placeholder_endpoint(atoms) -> bool:
    """Return whether an endpoint carries ATST placeholder results."""
    return atoms.info.get(ENDPOINT_RESULT_KEY) == ENDPOINT_PLACEHOLDER


def get_endpoint_results(atoms) -> tuple[float, np.ndarray] | None:
    """Return existing endpoint energy and forces without starting a calculation."""
    if is_placeholder_endpoint(atoms):
        return None
    # ASE owns validity of cached E/F for the current atoms. The CP geometry
    # fingerprint is a restart checkpoint guard, not an NEB admission proof.
    if CP_FACTS_INFO_KEY in atoms.info and read_constant_potential_facts(atoms, check_geometry=False) is None:
        return None
    calculator = atoms.calc
    results = getattr(calculator, "results", None)
    if not isinstance(results, dict) or not {"energy", "forces"}.issubset(results):
        return None
    try:
        if calculator.calculation_required(atoms, ["energy", "forces"]):
            return None
        energy = float(atoms.get_potential_energy())
        forces = np.asarray(atoms.get_forces(), dtype=float)
    except Exception:
        return None
    if not np.isfinite(energy) or forces.shape != (len(atoms), 3) or not np.all(np.isfinite(forces)):
        return None
    return energy, forces


def has_endpoint_results(atoms) -> bool:
    """Return whether an endpoint has usable non-placeholder results."""
    return get_endpoint_results(atoms) is not None


def has_trusted_endpoint_results(atoms) -> bool:
    """Compatibility alias for usable endpoint results regardless of origin."""
    return has_endpoint_results(atoms)


def expected_cp_endpoint_identity(config: dict[str, Any]) -> dict[str, Any] | None:
    """Return CP energy-definition context for NEB-family endpoint results."""
    cp_section = config.get("calculator", {}).get("constant_potential")
    if not isinstance(cp_section, dict):
        return None
    boundary = cp_section.get("energy_boundary", "reference_fcp")
    target = (
        cp_section.get("target_mu_ev")
        if boundary == "compensated_gate"
        else cp_section.get("potential_v")
    )
    return constant_potential_identity_for_config(config, target=target) if target is not None else None


def has_compatible_cp_endpoint_results(atoms, expected: dict[str, Any] | None) -> bool:
    """Check that endpoint E/F use the requested CP energy definition.

    CP endpoint energies need a complete CP result envelope and the same energy
    quantity as the requested run. Ordinary NEB cannot reuse a CP grand
    potential. Calculator settings and source provenance do not decide whether
    otherwise usable endpoint results may be reused.
    """
    endpoint_results = get_endpoint_results(atoms)
    if endpoint_results is None:
        return False
    if expected is None:
        return CP_FACTS_INFO_KEY not in atoms.info
    payload = read_constant_potential_facts(atoms, check_geometry=False)
    if payload is None:
        return False
    identity = payload.get("identity")
    facts = payload.get("facts")
    if not isinstance(facts, dict):
        return False
    try:
        facts_energy = float(facts["energy"])
    except (KeyError, TypeError, ValueError):
        return False
    if not np.isclose(facts_energy, endpoint_results[0], rtol=1e-10, atol=1e-12):
        return False
    return isinstance(identity, dict) and all(
        identity.get(key) == expected.get(key)
        for key in ("energy_boundary", "energy_definition")
    )


def endpoint_constant_potential_identity(atoms) -> dict[str, Any] | None:
    """Return the CP identity staged with an endpoint, when present."""
    facts = read_constant_potential_facts(atoms, check_geometry=False)
    if facts is not None and isinstance(facts.get("identity"), dict):
        return dict(facts["identity"])
    value = atoms.info.get(CP_ENDPOINT_IDENTITY_KEY)
    return dict(value) if isinstance(value, dict) else None


def constant_potential_identity_matches(atoms, expected: dict[str, Any] | None) -> bool:
    """Return whether staged CP facts carry the requested full identity.

    This legacy checkpoint helper is not used for NEB endpoint admission.
    """
    if expected is None:
        return True
    facts = read_constant_potential_facts(atoms)
    if facts is None or not isinstance(facts.get("identity"), dict):
        return False
    return _mapping_contains(facts["identity"], expected)


def freeze_current_results(atoms, status: str = ENDPOINT_COMPUTED):
    """Freeze current calculator results on an endpoint using SinglePointCalculator."""
    source_calculator = atoms.calc
    identity = getattr(source_calculator, "constant_potential_identity", None)
    energy = float(atoms.get_potential_energy())
    forces = np.asarray(atoms.get_forces(), dtype=float)
    kwargs: dict[str, Any] = {"energy": energy, "forces": forces}
    results = getattr(atoms.calc, "results", {}) if atoms.calc is not None else {}
    if "stress" in results:
        kwargs["stress"] = np.asarray(results["stress"], dtype=float)
    if identity is not None:
        # Explicitly copy the CP envelope before replacing the live calculator;
        # SinglePointCalculator does not retain backend-specific facts.
        if not isinstance(results, dict) or results.get("cp_converged") is not True:
            clear_constant_potential_facts(atoms)
        else:
            publish_constant_potential_facts(atoms, results, identity)
    else:
        clear_constant_potential_facts(atoms)
    atoms.calc = SinglePointCalculator(atoms, **kwargs)
    mark_endpoint_result(atoms, status)
    return atoms


def freeze_results(atoms, energy: float, forces, status: str | None, *, cp_facts: dict[str, Any] | None = None):
    """Freeze explicit endpoint results, preserving supplied CP facts when present."""
    clear_constant_potential_facts(atoms)
    atoms.calc = SinglePointCalculator(atoms, energy=float(energy), forces=np.asarray(forces, dtype=float))
    mark_endpoint_result(atoms, status)
    if cp_facts is not None:
        publish_constant_potential_facts(atoms, cp_facts["facts"], cp_facts["identity"])
    return atoms


def endpoint_policy(config: dict[str, Any], default: str = "auto") -> str:
    """Return and validate endpoint single-point policy from workflow config."""
    policy = config.get("endpoint_singlepoint", default)
    if policy not in {"auto", "always", "never"}:
        raise ValueError("endpoint_singlepoint must be one of: auto, always, never")
    return policy


def ensure_neb_endpoint_results(
    images: Sequence,
    get_calculator: Callable[[str], Any],
    policy: str = "auto",
    directories: tuple[str, str] = ("endpoint_initial", "endpoint_final"),
    context: str = "NEB",
    constant_potential_identity: dict[str, Any] | None = None,
):
    """Ensure the first and last images carry meaningful endpoint results.

    Args:
        images: NEB chain including endpoints.
        get_calculator: Callable receiving a directory suffix and returning an
            ASE calculator.
        policy: ``auto`` reuses readable finite endpoint energy/forces regardless
            of their origin and computes missing results with the current run's
            calculator; ``always`` recomputes both endpoints; ``never`` requires
            existing meaningful energy/force results.
        directories: Directory suffixes for initial and final endpoint
            calculations.
        context: Text used in warning/error messages.

    Returns:
        The input images sequence after in-place endpoint preparation.
    """
    if policy not in {"auto", "always", "never"}:
        raise ValueError("endpoint_singlepoint must be one of: auto, always, never")
    if len(images) < 2:
        raise ValueError("NEB endpoint preparation requires at least two images")

    endpoint_specs = ((0, "initial", directories[0]), (-1, "final", directories[1]))
    for index, label, directory in endpoint_specs:
        atoms = images[index]
        if policy == "never":
            if not has_endpoint_results(atoms):
                raise ValueError(
                    f"{context} {label} endpoint lacks meaningful energy/force results. "
                    "Relax the endpoint to obtain energy/forces, or allow endpoint "
                    "optimization before NEB."
                )
            if not has_compatible_cp_endpoint_results(atoms, constant_potential_identity):
                if constant_potential_identity is None:
                    raise ValueError(
                        f"{context} {label} endpoint carries a constant-potential "
                        "energy definition; set endpoint_singlepoint=auto to compute "
                        "ordinary endpoint energy/forces."
                    )
                raise ValueError(
                    f"{context} {label} endpoint lacks complete constant-potential "
                    "results with the requested energy definition; set "
                    "endpoint_singlepoint=auto to compute them."
                )
            continue  # never = explicit preserve of user-provided readable endpoints
        if policy == "auto" and has_compatible_cp_endpoint_results(atoms, constant_potential_identity):
            continue
        if policy == "auto":
            print(
                f"Warning: {context} {label} endpoint lacks usable energy/force "
                "results; running an endpoint single-point calculation."
            )
        elif policy == "always":
            print(f"Warning: {context} {label} endpoint is being recomputed by endpoint_singlepoint=always.")

        atoms.calc = get_calculator(directory)
        freeze_current_results(atoms, status=ENDPOINT_COMPUTED)
        if not has_compatible_cp_endpoint_results(atoms, constant_potential_identity):
            raise ValueError(f"{context} {label} endpoint calculation did not provide usable results")
    return images
