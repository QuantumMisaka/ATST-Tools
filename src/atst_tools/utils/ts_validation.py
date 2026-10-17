"""Evidence-limited transition-state screening from stored vibration results."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from numbers import Integral
from typing import Any

import numpy as np


SCHEMA_VERSION = "atst-ts-validation-v1"


def frozen_structure_fact(atoms: Any) -> dict[str, Any]:
    """Capture geometry, constraint scope and only matching frozen forces.

    Args:
        atoms: ASE ``Atoms`` object. Its calculator is never evaluated.

    Returns:
        JSON-ready evidence. ``forces`` is present only for a matching
        ``SinglePointCalculator`` snapshot. Unsupported constraint types leave
        ``fixed_indices`` unknown.
    """
    from ase.calculators.singlepoint import SinglePointCalculator
    from ase.constraints import FixAtoms

    fact = {
        "n_atoms": len(atoms),
        "geometry": {
            "numbers": atoms.numbers.tolist(),
            "positions": atoms.positions.tolist(),
            "cell": atoms.cell.array.tolist(),
            "pbc": atoms.pbc.tolist(),
        },
    }
    fixed = []
    for constraint in atoms.constraints:
        if not isinstance(constraint, FixAtoms):
            fixed = None
            break
        fixed.extend(int(index) for index in constraint.get_indices())
    fact["fixed_indices"] = sorted(set(fixed)) if fixed is not None else None
    calculator = atoms.calc
    if isinstance(calculator, SinglePointCalculator) and not calculator.check_state(atoms):
        forces = calculator.results.get("forces")
        if forces is not None:
            fact["forces"] = np.asarray(forces).tolist()
    return fact


def _finite_number(value: Any, name: str) -> float:
    """Convert one numeric fact to a finite float with a useful error."""
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(f"{name} must be a finite number, got {value!r}.")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"{name} must be a finite number, got {value!r}.") from exc
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite, got {value!r}.")
    return number


def _indices(value: Any, name: str, n_atoms: int | None, *, allow_empty: bool = False) -> list[int] | None:
    """Validate atom indices without interpreting text or booleans as indices."""
    if value is None:
        return None
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise TypeError(f"{name} must be a sequence of atom indices.")
    result = []
    for item in value:
        if isinstance(item, (bool, np.bool_)) or not isinstance(item, Integral):
            raise TypeError(f"{name} entries must be integer atom indices, got {item!r}.")
        index = int(item)
        if index < 0 or (n_atoms is not None and index >= n_atoms):
            raise ValueError(f"{name} index {index} is out of range.")
        if index in result:
            raise ValueError(f"{name} contains duplicate index {index}.")
        result.append(index)
    if not result and not allow_empty:
        raise ValueError(f"{name} must not be empty.")
    return result


def _atom_count(value: Any) -> int | None:
    """Validate an optional atom count."""
    if value is None:
        return None
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral):
        raise TypeError("n_atoms must be a positive integer.")
    if value <= 0:
        raise ValueError("n_atoms must be a positive integer.")
    return int(value)


def _stored_forces(value: Any, n_atoms: int | None) -> list[list[float]] | None:
    """Validate frozen Cartesian forces without calling an ASE calculator."""
    if value is None:
        return None
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise TypeError("forces must be a sequence of three-component vectors.")
    if n_atoms is None or len(value) != n_atoms:
        raise ValueError("forces must have one three-component vector per atom.")
    forces = []
    for index, vector in enumerate(value):
        if isinstance(vector, (str, bytes)) or not isinstance(vector, Sequence) or len(vector) != 3:
            raise ValueError(f"forces[{index}] must have three components.")
        forces.append([_finite_number(item, f"forces[{index}]") for item in vector])
    return forces


def _geometry(value: Any, n_atoms: int | None) -> dict[str, Any] | None:
    """Validate a complete frozen geometry before pairing force evidence."""
    if value is None:
        return None
    if not isinstance(value, Mapping):
        raise TypeError("geometry must be a mapping.")
    if n_atoms is None:
        raise ValueError("geometry requires n_atoms.")
    numbers = value.get("numbers")
    positions = value.get("positions")
    cell = value.get("cell")
    pbc = value.get("pbc")
    if not isinstance(numbers, Sequence) or isinstance(numbers, (str, bytes)) or len(numbers) != n_atoms:
        raise ValueError("geometry.numbers must have one atomic number per atom.")
    if any(isinstance(number, (bool, np.bool_)) or not isinstance(number, Integral) or number <= 0 for number in numbers):
        raise ValueError("geometry.numbers must contain positive integers.")
    if not isinstance(positions, Sequence) or isinstance(positions, (str, bytes)) or len(positions) != n_atoms:
        raise ValueError("geometry.positions must have one vector per atom.")
    if not isinstance(cell, Sequence) or isinstance(cell, (str, bytes)) or len(cell) != 3:
        raise ValueError("geometry.cell must have three vectors.")
    if not isinstance(pbc, Sequence) or isinstance(pbc, (str, bytes)) or len(pbc) != 3:
        raise ValueError("geometry.pbc must have three booleans.")
    if any(not isinstance(flag, (bool, np.bool_)) for flag in pbc):
        raise TypeError("geometry.pbc must contain booleans.")
    return {
        "numbers": [int(number) for number in numbers],
        "positions": _stored_forces(positions, n_atoms),
        "cell": _stored_forces(cell, 3),
        "pbc": [bool(flag) for flag in pbc],
    }


def build_ts_validation_summary(
    vibration_results: Mapping[str, Any],
    *,
    fmax: float | None = None,
    fmax_threshold: float = 0.05,
    imaginary_cutoff_cm1: float = 50.0,
    source: str | None = None,
    structure: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Screen stored vibration and force facts on their declared atom scope.

    The historical ``status`` and ``checks`` fields remain a vibration/optional
    scalar-fmax screen. ``requested_checks_status`` includes evidence coverage
    and stationarity, and never asserts global saddle or IRC connectivity.

    Args:
        vibration_results: Stored vibration output, including optional atom indices.
        fmax: Legacy caller-supplied maximum force for the historical check.
        fmax_threshold: Maximum force threshold in eV/Angstrom.
        imaginary_cutoff_cm1: Significant imaginary frequency cutoff in cm^-1.
        source: Vibration result source identifier.
        structure: Frozen structure facts with ``n_atoms``, ``fixed_indices``
            and optional stored ``forces``. No calculator is evaluated.

    Returns:
        JSON-ready scope-limited screening summary.

    Raises:
        TypeError: A supplied fact has the wrong type.
        ValueError: A supplied numeric fact is nonfinite or structurally invalid.
    """
    if not isinstance(vibration_results, Mapping):
        raise TypeError("vibration_results must be a mapping.")
    if structure is not None and not isinstance(structure, Mapping):
        raise TypeError("structure must be a mapping of frozen facts.")
    structure = structure or {}
    cutoff = _finite_number(imaginary_cutoff_cm1, "imaginary_cutoff_cm1")
    threshold = _finite_number(fmax_threshold, "fmax_threshold")
    if cutoff < 0 or threshold < 0:
        raise ValueError("Cutoffs and thresholds must be non-negative.")

    raw_imaginary = vibration_results.get("imaginary_frequencies")
    if raw_imaginary is None:
        imag_values = None
    elif isinstance(raw_imaginary, (str, bytes)) or not isinstance(raw_imaginary, Sequence):
        raise TypeError("imaginary_frequencies must be a sequence of finite numbers.")
    else:
        imag_values = [abs(_finite_number(value, "imaginary_frequencies entry")) for value in raw_imaginary]
    significant = None if imag_values is None else [value for value in imag_values if value >= cutoff]
    n_imag = None if significant is None else len(significant)
    imag_ok = None if n_imag is None else n_imag == 1

    result_count = _atom_count(vibration_results.get("n_atoms"))
    structure_count = _atom_count(structure.get("n_atoms"))
    if result_count is not None and structure_count is not None and result_count != structure_count:
        raise ValueError("Vibration and structure n_atoms disagree.")
    n_atoms = result_count if result_count is not None else structure_count
    result_geometry = _geometry(vibration_results.get("geometry"), n_atoms)
    structure_geometry = _geometry(structure.get("geometry"), n_atoms)
    if result_geometry is not None and structure_geometry is not None and result_geometry != structure_geometry:
        raise ValueError("Frozen structure geometry does not match vibration geometry.")
    geometry_paired = result_geometry is not None and structure_geometry is not None
    selected = _indices(vibration_results.get("indices"), "indices", n_atoms)
    fixed_source = vibration_results.get("fixed_indices")
    if geometry_paired and fixed_source is not None and structure.get("fixed_indices") is not None:
        if set(_indices(fixed_source, "fixed_indices", n_atoms, allow_empty=True)) != set(
            _indices(structure["fixed_indices"], "fixed_indices", n_atoms, allow_empty=True)
        ):
            raise ValueError("Frozen structure constraints disagree with vibration constraints.")
    if fixed_source is None and geometry_paired:
        fixed_source = structure.get("fixed_indices")
    fixed = _indices(fixed_source, "fixed_indices", n_atoms, allow_empty=True)
    fixed_in_selection = (
        None if selected is None or fixed is None
        else [index for index in selected if index in set(fixed)]
    )
    # A missing index list does not establish which atoms ASE vibrated. Older
    # result files may omit this fact even when atom count and geometry exist.
    active = selected
    free = None if n_atoms is None or fixed is None else [i for i in range(n_atoms) if i not in fixed]
    covers_free = None if free is None or active is None else set(free).issubset(active)
    full_system = n_atoms is not None and fixed == [] and active is not None and set(active) == set(range(n_atoms))
    if full_system:
        kind = "full_system"
    elif fixed and selected is not None:
        kind = "constrained"
    elif selected is not None and n_atoms is not None and len(selected) < n_atoms:
        kind = "subset"
    else:
        kind = "unknown"

    # A force array without matching geometry is not evidence for the
    # vibrated structure. Never read calculator state or launch an evaluation.
    supplied_forces = _stored_forces(structure.get("forces"), n_atoms)
    forces = supplied_forces if geometry_paired else None
    # ASE permits explicit selection of fixed atoms. Their stored Cartesian
    # forces do not establish stationarity for that mixed scope.
    if forces is not None and active is not None and not fixed_in_selection:
        measured_fmax = max((math.hypot(*forces[i]) for i in active), default=None)
    else:
        measured_fmax = None
    stationary = None if measured_fmax is None else measured_fmax <= threshold
    stationarity_status = "unknown" if stationary is None else ("pass" if stationary else "fail")

    legacy_fmax = None if fmax is None else _finite_number(fmax, "fmax")
    legacy_fmax_ok = None if legacy_fmax is None else legacy_fmax <= threshold
    legacy_status = "pass" if imag_ok and legacy_fmax_ok is not False else "fail"
    if imag_ok is False or stationary is False:
        requested_status = "fail"
    elif imag_ok is None or stationary is None or covers_free is not True or fixed_in_selection:
        requested_status = "unknown"
    else:
        requested_status = "pass"

    return {
        "schema_version": SCHEMA_VERSION,
        "source": source,
        "status": legacy_status,
        "checks": {
            "n_imaginary_modes": {"value": n_imag, "expected": 1, "imaginary_cutoff_cm1": imaginary_cutoff_cm1, "pass": imag_ok},
            "fmax": {"value": fmax, "threshold": fmax_threshold, "pass": legacy_fmax_ok},
        },
        "imaginary_frequencies": significant if significant is not None else [],
        "scope": {
            "n_atoms": n_atoms,
            "indices": selected,
            "fixed_indices": fixed,
            "fixed_indices_in_selection": fixed_in_selection,
            "active_indices": active,
            "kind": kind,
            "full_system": full_system,
            "covers_all_unconstrained_atoms": covers_free,
            "geometry_paired": geometry_paired,
        },
        "stationarity": {"status": stationarity_status, "fmax": measured_fmax, "threshold": threshold, "scope": "active_indices"},
        "requested_checks_status": requested_status,
    }
