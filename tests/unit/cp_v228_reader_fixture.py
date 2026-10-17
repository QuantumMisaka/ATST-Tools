"""Frozen v2.2.8 constant-potential facts reader used for compatibility tests.

The definitions below are copied from tag v2.2.8 at
src/atst_tools/calculators/constant_potential.py. This fixture contains only
the reader, restart validation, and their direct dependencies. Do not update
it to match current production behavior.
"""

from __future__ import annotations

from collections.abc import Mapping
import copy
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

_ENERGY_BOUNDARIES = {"reference_fcp", "compensated_gate"}


class ConstantPotentialError(RuntimeError):
    """Base exception for constant-potential evaluation failures."""


class ConstantPotentialEvaluationError(ConstantPotentialError):
    """Raised when one fixed-electron backend evaluation is unusable."""


CP_FACTS_INFO_KEY = "atst_constant_potential_facts"


CP_IDENTITY_INFO_KEY = "atst_constant_potential_identity"


CP_FACTS_SCHEMA = "atst-constant-potential-facts-v1"


_CP_FACT_KEYS = (
    "energy",
    "free_energy",
    "raw_energy",
    "raw_free_energy",
    "efermi",
    "vacuum_level",
    "fermishift",
    "energy_boundary",
    "profile_status",
    "boundary_parameters",
    "compensation",
    "mu_calc",
    "mu_target",
    "potential_calc",
    "potential_target",
    "residual_mu",
    "residual_v",
    "omega_correction",
    "nelec",
    "reference_electrons",
    "delta_nelec",
    "forces",
    "scf_converged",
    "cp_converged",
    "cp_iterations",
    "cp_history",
)


_CP_REQUIRED_FACT_KEYS = frozenset(_CP_FACT_KEYS)


_CP_REQUIRED_IDENTITY_KEYS = frozenset(
    {
        "energy_boundary",
        "potential_v",
        "target_mu_ev",
        "reference_electrode",
        "work_ref",
        "work_ref_source",
        "reference_pH",
        "temperature_K",
        "reference_electrons",
        "boundary_parameters",
        "energy_definition",
        "calculator",
    }
)


_CP_REQUIRED_FINITE_FACT_KEYS = frozenset(
    {
        "energy",
        "free_energy",
        "raw_energy",
        "raw_free_energy",
        "efermi",
        "mu_calc",
        "mu_target",
        "residual_mu",
        "residual_v",
        "omega_correction",
        "nelec",
        "reference_electrons",
        "delta_nelec",
    }
)


def _jsonable(value: Any) -> Any:
    """Convert common ASE/NumPy values into JSON-safe values."""
    if isinstance(value, Mapping):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.integer, np.floating, np.bool_)):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    return value


def constant_potential_geometry_fingerprint(atoms: Any) -> str:
    """Return a deterministic fingerprint of the ASE geometry and cell."""
    digest = hashlib.sha256()
    numbers = np.asarray(atoms.numbers, dtype="<i8")
    positions = np.asarray(atoms.positions, dtype="<f8")
    cell = np.asarray(atoms.cell.array, dtype="<f8")
    pbc = np.asarray(atoms.pbc, dtype=np.uint8)
    for value in (numbers, positions, cell, pbc):
        digest.update(np.asarray(value).tobytes(order="C"))
    return digest.hexdigest()


def _mapping_contains(actual: Any, expected: Any) -> bool:
    """Compare an expected identity as a recursive subset of actual facts."""
    if isinstance(expected, Mapping):
        if not isinstance(actual, Mapping):
            return False
        return all(
            key == "provenance"
            or (key in actual and _mapping_contains(actual[key], value))
            for key, value in expected.items()
        )
    if isinstance(expected, (list, tuple)):
        if not isinstance(actual, (list, tuple)) or len(actual) != len(expected):
            return False
        return all(_mapping_contains(item, value) for item, value in zip(actual, expected))
    return actual == expected


def _same_finite_float(left: Any, right: Any) -> bool:
    """Compare two serialized scalar values without accepting non-finite data."""
    try:
        left_value = float(left)
        right_value = float(right)
    except (TypeError, ValueError):
        return False
    if not np.isfinite(left_value) or not np.isfinite(right_value):
        return False
    return bool(np.isclose(left_value, right_value, rtol=1e-10, atol=1e-12))


def _constant_potential_facts_error(
    atoms: Any, payload: Any, *, require_duplicate_identity: bool = True
) -> str | None:
    """Return a fail-closed validation error for one durable CP envelope."""
    if atoms is None or not hasattr(atoms, "info"):
        return "constant-potential facts require an ASE atoms object"
    if not isinstance(payload, Mapping):
        return "constant-potential facts envelope is not a mapping"
    if payload.get("schema") != CP_FACTS_SCHEMA or payload.get("workflow") != "constant_potential":
        return "constant-potential facts schema or workflow is unsupported"
    if payload.get("status") != "complete" or payload.get("cp_converged") is not True:
        return "constant-potential facts are not a complete converged result"
    try:
        geometry_fingerprint = constant_potential_geometry_fingerprint(atoms)
    except (AttributeError, TypeError, ValueError):
        return "constant-potential facts geometry cannot be fingerprinted"
    if payload.get("geometry_fingerprint") != geometry_fingerprint:
        return "constant-potential facts geometry fingerprint does not match"

    identity = payload.get("identity")
    facts = payload.get("facts")
    if not isinstance(identity, Mapping):
        return "constant-potential facts identity is missing or malformed"
    if not isinstance(facts, Mapping):
        return "constant-potential facts facts mapping is missing or malformed"
    duplicate_identity = atoms.info.get(CP_IDENTITY_INFO_KEY)
    if require_duplicate_identity and (
        not isinstance(duplicate_identity, Mapping)
        or _jsonable(duplicate_identity) != _jsonable(identity)
    ):
        return "constant-potential facts duplicate identity is missing or inconsistent"

    missing_identity = sorted(key for key in _CP_REQUIRED_IDENTITY_KEYS if key not in identity)
    if missing_identity:
        return f"constant-potential facts identity is missing required fields: {', '.join(missing_identity)}"
    missing_facts = sorted(key for key in _CP_REQUIRED_FACT_KEYS if key not in facts)
    if missing_facts:
        return f"constant-potential facts are missing required fields: {', '.join(missing_facts)}"

    if facts["cp_converged"] is not True or facts["scf_converged"] is not True:
        return "constant-potential facts require cp_converged=true and scf_converged=true"
    if not isinstance(facts["energy_boundary"], str) or facts["energy_boundary"] not in _ENERGY_BOUNDARIES:
        return "constant-potential facts energy_boundary is unsupported"
    if facts["energy_boundary"] != identity["energy_boundary"]:
        return "constant-potential facts boundary disagrees with identity"
    if not isinstance(facts["profile_status"], str) or not facts["profile_status"]:
        return "constant-potential facts profile_status is malformed"
    if not isinstance(facts["boundary_parameters"], Mapping):
        return "constant-potential facts boundary_parameters is malformed"
    if not isinstance(identity["boundary_parameters"], Mapping):
        return "constant-potential facts identity boundary_parameters is malformed"
    if not _mapping_contains(identity["boundary_parameters"], facts["boundary_parameters"]):
        return "constant-potential facts boundary parameters disagree with identity"
    if not isinstance(identity["calculator"], str) or not identity["calculator"]:
        return "constant-potential facts identity calculator is malformed"
    expected_definition = (
        "omega_compensated_gate"
        if facts["energy_boundary"] == "compensated_gate"
        else "omega_aligned"
    )
    if identity["energy_definition"] != expected_definition:
        return "constant-potential facts identity energy definition is inconsistent"

    for key in _CP_REQUIRED_FINITE_FACT_KEYS:
        if not _same_finite_float(facts[key], facts[key]):
            return f"constant-potential facts {key} is not finite"
    if not _same_finite_float(identity["reference_electrons"], identity["reference_electrons"]):
        return "constant-potential facts identity reference_electrons is not finite"
    if not _same_finite_float(facts["reference_electrons"], identity["reference_electrons"]):
        return "constant-potential facts reference_electrons disagrees with identity"
    if not _same_finite_float(
        facts["delta_nelec"], float(facts["nelec"]) - float(facts["reference_electrons"])
    ):
        return "constant-potential facts delta_nelec is inconsistent with electron counts"
    if not _same_finite_float(
        facts["residual_mu"], float(facts["mu_calc"]) - float(facts["mu_target"])
    ):
        return "constant-potential facts residual_mu is inconsistent with mu_calc and mu_target"
    if not _same_finite_float(
        facts["energy"], float(facts["raw_energy"]) + float(facts["omega_correction"])
    ) or not _same_finite_float(
        facts["free_energy"], float(facts["raw_free_energy"]) + float(facts["omega_correction"])
    ):
        return "constant-potential facts energy correction is inconsistent"

    if facts["energy_boundary"] == "compensated_gate":
        if not _same_finite_float(identity["target_mu_ev"], facts["mu_target"]):
            return "constant-potential facts mu_target disagrees with identity"
        if facts["vacuum_level"] is not None or facts["fermishift"] is not None:
            return "constant-potential compensated-gate facts must not contain a vacuum reference"
        if facts["potential_calc"] is not None or facts["potential_target"] is not None:
            return "constant-potential compensated-gate facts must not contain a voltage target"
        if not _same_finite_float(facts["residual_v"], facts["residual_mu"]):
            return "constant-potential compensated-gate residual_v must equal residual_mu"
        compensation = facts["compensation"]
        if not isinstance(compensation, Mapping):
            return "constant-potential compensated-gate compensation evidence is missing"
        required_compensation = (
            "boundary",
            "gate_derivative_ev",
            "dipole_derivative_ev",
            "compensation_derivative_ev",
            "integrated_electrons",
            "density_electron_error",
            "ionic_valence_electrons",
            "total_dipole_ry_au",
            "density_precision",
            "density_electron_tolerance",
        )
        missing_compensation = [
            key for key in required_compensation if key not in compensation
        ]
        if missing_compensation:
            return (
                "constant-potential compensated-gate compensation is missing required fields: "
                + ", ".join(missing_compensation)
            )
        if compensation["boundary"] != "compensated_gate":
            return "constant-potential compensation boundary is inconsistent"
        for key in (
            "gate_derivative_ev",
            "dipole_derivative_ev",
            "compensation_derivative_ev",
            "integrated_electrons",
            "density_electron_error",
            "ionic_valence_electrons",
            "total_dipole_ry_au",
            "density_electron_tolerance",
        ):
            if not _same_finite_float(compensation[key], compensation[key]):
                return f"constant-potential compensation {key} is not finite"
        try:
            density_precision = float(compensation["density_precision"])
        except (TypeError, ValueError):
            return "constant-potential compensation density_precision is malformed"
        if not np.isfinite(density_precision) or not density_precision.is_integer() or density_precision < 10:
            return "constant-potential compensation density_precision is invalid"
        if float(compensation["density_electron_tolerance"]) <= 0:
            return "constant-potential compensation density_electron_tolerance is invalid"
        if not _same_finite_float(
            compensation["compensation_derivative_ev"],
            float(compensation["gate_derivative_ev"])
            + float(compensation["dipole_derivative_ev"]),
        ):
            return "constant-potential compensation derivative is inconsistent"
        if not _same_finite_float(
            facts["mu_calc"],
            float(facts["efermi"]) + float(compensation["compensation_derivative_ev"]),
        ):
            return "constant-potential compensated chemical potential is inconsistent"
        if not _same_finite_float(
            facts["omega_correction"], -float(facts["mu_target"]) * float(facts["delta_nelec"])
        ):
            return "constant-potential compensated grand-potential correction is inconsistent"
        # Match the same PP-valence tolerance used by read_gate_compensation.
        if abs(float(compensation["ionic_valence_electrons"]) - float(facts["reference_electrons"])) > 1e-6:
            return "constant-potential compensation ionic valence disagrees with reference electrons"
        if not _same_finite_float(compensation["density_electron_error"],
                                  float(compensation["integrated_electrons"])
                                  - float(facts["nelec"])):
            return "constant-potential compensation density error is inconsistent"
        if abs(float(compensation["density_electron_error"])) > float(
            compensation["density_electron_tolerance"]
        ):
            return "constant-potential compensation density error exceeds tolerance"
    else:
        for key in ("vacuum_level", "fermishift", "potential_calc", "potential_target"):
            if not _same_finite_float(facts[key], facts[key]):
                return f"constant-potential reference fact {key} is not finite"
        if not _same_finite_float(facts["fermishift"], -float(facts["vacuum_level"])):
            return "constant-potential reference fermishift is inconsistent with vacuum_level"
        if not _same_finite_float(
            facts["mu_calc"], float(facts["efermi"]) + float(facts["fermishift"])
        ):
            return "constant-potential reference chemical potential is inconsistent"
        if not _same_finite_float(
            facts["omega_correction"],
            (float(facts["fermishift"]) - float(facts["mu_target"])) * float(facts["delta_nelec"]),
        ):
            return "constant-potential reference grand-potential correction is inconsistent"
        if not _same_finite_float(identity["potential_v"], facts["potential_target"]):
            return "constant-potential facts potential_target disagrees with identity"
        try:
            expected_mu = -(
                float(identity["potential_v"]) + float(identity["work_ref"])
            )
        except (TypeError, ValueError):
            return "constant-potential reference identity work_ref is malformed"
        if not _same_finite_float(facts["mu_target"], expected_mu):
            return "constant-potential facts mu_target disagrees with identity"
        try:
            expected_residual_v = float(facts["potential_calc"]) - float(facts["potential_target"])
        except (TypeError, ValueError):
            return "constant-potential reference potential values are malformed"
        if not _same_finite_float(facts["residual_v"], expected_residual_v):
            return "constant-potential facts residual_v is inconsistent with voltage values"
    try:
        forces = np.asarray(facts["forces"], dtype=float)
    except (TypeError, ValueError):
        return "constant-potential facts forces are malformed"
    if forces.shape != (len(atoms), 3) or not np.all(np.isfinite(forces)):
        return "constant-potential facts forces have invalid shape or values"
    try:
        iterations = float(facts["cp_iterations"])
    except (TypeError, ValueError):
        return "constant-potential facts cp_iterations is malformed"
    if not np.isfinite(iterations) or iterations < 0 or not iterations.is_integer():
        return "constant-potential facts cp_iterations is invalid"
    if not isinstance(facts["cp_history"], list) or not facts["cp_history"]:
        return "constant-potential facts cp_history is missing or empty"
    try:
        json.dumps(_jsonable(payload), allow_nan=False)
    except (TypeError, ValueError):
        return "constant-potential facts are not finite JSON"
    return None


def publish_constant_potential_facts(
    atoms: Any, results: Mapping[str, Any], identity: Mapping[str, Any]
) -> dict[str, Any]:
    """Persist one validated CP result envelope on the evaluated atoms."""
    if not isinstance(results, Mapping) or not isinstance(identity, Mapping):
        raise ConstantPotentialEvaluationError("constant-potential facts require result and identity mappings")
    facts = {key: _jsonable(results[key]) for key in _CP_FACT_KEYS if key in results}
    payload = {
        "schema": CP_FACTS_SCHEMA,
        "workflow": "constant_potential",
        "status": "complete",
        "geometry_fingerprint": constant_potential_geometry_fingerprint(atoms),
        "identity": _jsonable(identity),
        "facts": facts,
        "cp_converged": facts.get("cp_converged") is True,
    }
    try:
        json.dumps(payload, allow_nan=False)
    except (TypeError, ValueError) as exc:
        raise ConstantPotentialEvaluationError("constant-potential facts are not finite JSON") from exc
    validation_error = _constant_potential_facts_error(
        atoms, payload, require_duplicate_identity=False
    )
    if validation_error is not None:
        raise ConstantPotentialEvaluationError(validation_error)
    atoms.info[CP_FACTS_INFO_KEY] = payload
    atoms.info[CP_IDENTITY_INFO_KEY] = copy.deepcopy(payload["identity"])
    return copy.deepcopy(payload)


def read_constant_potential_facts(atoms: Any) -> dict[str, Any] | None:
    """Read and validate CP facts persisted on an ASE atoms object."""
    if atoms is None or not hasattr(atoms, "info"):
        return None
    payload = atoms.info.get(CP_FACTS_INFO_KEY)
    if _constant_potential_facts_error(atoms, payload) is not None:
        return None
    return copy.deepcopy(dict(payload))


def constant_potential_restart_initial_electrons(
    atoms: Any,
    expected_identity: Mapping[str, Any],
    potential_tolerance_v: float,
    *,
    allow_recompute: bool = False,
) -> float | None:
    """Return the last evaluated electron count for a matching CP frame.

    A structure trajectory is a nuclear checkpoint; it does not serialize the
    calculator's custom result fields.  This helper consumes the durable facts
    envelope copied into ``Atoms.info`` and therefore only returns an electron
    count after validating the complete SCF/CP result, geometry, fixed
    Hamiltonian identity, target, and current residual tolerance.  A missing or
    stale envelope may be treated as an explicit recomputation when
    ``allow_recompute`` is true; otherwise restart fails closed.

    The returned value is the completed evaluation's ``nelec``.  It never uses
    a Newton candidate or a partially written history record.
    """
    if not isinstance(expected_identity, Mapping):
        raise ConstantPotentialError("constant-potential restart identity is missing or malformed")
    try:
        tolerance = float(potential_tolerance_v)
    except (TypeError, ValueError) as exc:
        raise ConstantPotentialError("constant-potential restart tolerance is malformed") from exc
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ConstantPotentialError("constant-potential restart tolerance must be positive and finite")

    has_payload = atoms is not None and hasattr(atoms, "info") and CP_FACTS_INFO_KEY in atoms.info
    if not has_payload:
        if allow_recompute:
            return None
        raise ConstantPotentialError("constant-potential restart checkpoint is missing")
    raw_payload = atoms.info.get(CP_FACTS_INFO_KEY)
    validation_error = _constant_potential_facts_error(atoms, raw_payload)
    payload = read_constant_potential_facts(atoms)
    if validation_error is not None or payload is None:
        if allow_recompute:
            return None
        detail = validation_error or "constant-potential facts could not be read"
        raise ConstantPotentialError(
            f"constant-potential restart checkpoint is stale or malformed: {detail}"
        )

    identity = payload.get("identity")
    if not isinstance(identity, Mapping) or not _mapping_contains(identity, expected_identity):
        if allow_recompute:
            return None
        raise ConstantPotentialError(
            "constant-potential restart checkpoint identity does not match the requested calculation"
        )
    facts = payload.get("facts")
    if not isinstance(facts, Mapping):
        if allow_recompute:
            return None
        raise ConstantPotentialError("constant-potential restart checkpoint facts are malformed")

    residual_key = "residual_mu" if facts.get("energy_boundary") == "compensated_gate" else "residual_v"
    try:
        residual = float(facts[residual_key])
        nelec = float(facts["nelec"])
    except (KeyError, TypeError, ValueError) as exc:
        if allow_recompute:
            return None
        raise ConstantPotentialError(
            "constant-potential restart checkpoint is missing its evaluated residual or electron count"
        ) from exc
    if not np.isfinite(residual) or not np.isfinite(nelec) or nelec < 0:
        if allow_recompute:
            return None
        raise ConstantPotentialError(
            "constant-potential restart checkpoint has a non-finite or negative evaluated electron count"
        )
    if abs(residual) > tolerance:
        if allow_recompute:
            return None
        raise ConstantPotentialError(
            "constant-potential restart checkpoint residual exceeds the current tolerance"
        )
    return nelec
