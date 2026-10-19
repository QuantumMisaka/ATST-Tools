"""Compare constant-potential result records from two implementations.

The reference side may be produced by FCP v2 and ase-abacus; the candidate
side may be produced by ATST and abacuslite. The comparison is intentionally
field-oriented so a real-engine acceptance run can fail on one quantity
without losing the rest of the evidence.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys
from typing import Any, Mapping


DEFAULT_TOLERANCES: dict[str, float] = {
    "nelec": 1e-6,
    "efermi_eV": 1e-5,
    "vacuum_level_eV": 1e-5,
    "fermishift_eV": 1e-5,
    "mu_calc_eV": 1e-5,
    "mu_target_eV": 1e-6,
    "residual_mu_eV": 1e-5,
    "energy_eV": 1e-6,
    "free_energy_eV": 1e-6,
    "forces_eV_per_A": 1e-5,
    "density_integral_e": 1e-6,
    "compensation_derivative_eV": 1e-4,
    "capacitance_e_per_V": 1e-6,
}

REQUIRED_FIELDS = (
    "nelec",
    "efermi_eV",
    "mu_calc_eV",
    "energy_eV",
    "forces_eV_per_A",
)

INFORMATIONAL_FIELDS = frozenset({"iterations", "capacitance_e_per_V"})


def _flatten_numeric(value: Any, path: str, out: dict[str, Any]) -> None:
    """Flatten nested numeric sequences into path-keyed scalar entries."""
    if isinstance(value, (list, tuple)):
        for index, item in enumerate(value):
            _flatten_numeric(item, f"{path}[{index}]", out)
        return
    out[path] = value


def _mismatch(
    mismatches: list[dict[str, Any]],
    field: str,
    reference: Any,
    candidate: Any,
    atol: float | None = None,
    delta: float | None = None,
) -> None:
    """Append one explicit comparison failure."""
    mismatches.append(
        {
            "field": field,
            "reference": reference,
            "candidate": candidate,
            "atol": atol,
            "delta": delta,
        }
    )


def _compare_scalar(
    field: str,
    reference: Any,
    candidate: Any,
    atol: float,
    mismatches: list[dict[str, Any]],
) -> None:
    """Compare one optional scalar field."""
    if reference is None and candidate is None:
        return
    if reference is None or candidate is None:
        _mismatch(mismatches, field, reference, candidate, atol=atol)
        return
    try:
        reference_value = float(reference)
        candidate_value = float(candidate)
    except (TypeError, ValueError):
        _mismatch(mismatches, field, reference, candidate, atol=atol)
        return
    if not (
        math.isfinite(reference_value)
        and math.isfinite(candidate_value)
        and math.isclose(reference_value, candidate_value, rel_tol=0.0, abs_tol=atol)
    ):
        _mismatch(
            mismatches,
            field,
            reference_value,
            candidate_value,
            atol=atol,
            delta=abs(reference_value - candidate_value),
        )


def _compare_numeric_tree(
    field: str,
    reference: Any,
    candidate: Any,
    atol: float,
    mismatches: list[dict[str, Any]],
) -> None:
    """Compare a nested numeric sequence, including its shape."""
    reference_flat: dict[str, Any] = {}
    candidate_flat: dict[str, Any] = {}
    _flatten_numeric(reference, field, reference_flat)
    _flatten_numeric(candidate, field, candidate_flat)
    if set(reference_flat) != set(candidate_flat):
        _mismatch(
            mismatches,
            f"{field}.shape",
            sorted(reference_flat),
            sorted(candidate_flat),
            atol=atol,
        )
        return
    for path in sorted(reference_flat):
        _compare_scalar(path, reference_flat[path], candidate_flat[path], atol, mismatches)


def compare_records(
    reference: Mapping[str, Any],
    candidate: Mapping[str, Any],
    tolerances: Mapping[str, float] | None = None,
) -> dict[str, Any]:
    """Compare two constant-potential result mappings.

    Args:
        reference: Reference record, for example from FCP v2.
        candidate: Candidate record, for example from ATST.
        tolerances: Optional absolute tolerances overriding ``DEFAULT_TOLERANCES``.

    Returns:
        A JSON-compatible report with ``ok``, ``compared`` and ``mismatches``.

    Raises:
        ValueError: If a required field is absent from either record.
    """
    missing = [
        f"{side}:{field}"
        for side, record in (("reference", reference), ("candidate", candidate))
        for field in REQUIRED_FIELDS
        if field not in record
    ]
    if missing:
        raise ValueError("missing required fields: " + ", ".join(missing))

    applied = dict(DEFAULT_TOLERANCES)
    if tolerances:
        applied.update({str(key): float(value) for key, value in tolerances.items()})

    mismatches: list[dict[str, Any]] = []
    informational: list[dict[str, Any]] = []
    compared: list[str] = []
    for field in (
        "nelec",
        "efermi_eV",
        "vacuum_level_eV",
        "fermishift_eV",
        "mu_calc_eV",
        "mu_target_eV",
        "residual_mu_eV",
        "energy_eV",
        "free_energy_eV",
        "density_integral_e",
        "compensation_derivative_eV",
        "capacitance_e_per_V",
    ):
        if field not in reference and field not in candidate:
            continue
        compared.append(field)
        target = informational if field in INFORMATIONAL_FIELDS else mismatches
        _compare_scalar(
            field,
            reference.get(field),
            candidate.get(field),
            applied.get(field, 0.0),
            target,
        )

    compared.append("forces_eV_per_A")
    _compare_numeric_tree(
        "forces_eV_per_A",
        reference["forces_eV_per_A"],
        candidate["forces_eV_per_A"],
        applied["forces_eV_per_A"],
        mismatches,
    )

    for field in ("iterations", "converged"):
        if field in reference or field in candidate:
            compared.append(field)
            if reference.get(field) != candidate.get(field):
                _mismatch(
                    informational if field in INFORMATIONAL_FIELDS else mismatches,
                    field,
                    reference.get(field),
                    candidate.get(field),
                )

    return {
        "ok": not mismatches,
        "compared": compared,
        "tolerances": applied,
        "mismatches": mismatches,
        "informational": informational,
    }


def main(argv: list[str] | None = None) -> int:
    """Run the command-line comparison."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, help="Reference JSON record.")
    parser.add_argument("--candidate", required=True, help="Candidate JSON record.")
    parser.add_argument(
        "--tolerances",
        help="Optional JSON mapping of field names to absolute tolerances.",
    )
    args = parser.parse_args(argv)

    reference = json.loads(Path(args.reference).read_text(encoding="utf-8"))
    candidate = json.loads(Path(args.candidate).read_text(encoding="utf-8"))
    tolerances = (
        json.loads(Path(args.tolerances).read_text(encoding="utf-8"))
        if args.tolerances
        else None
    )
    try:
        report = compare_records(reference, candidate, tolerances)
    except ValueError as exc:
        print(json.dumps({"ok": False, "error": str(exc)}, indent=2))
        return 2
    print(json.dumps(report, indent=2))
    return 0 if report["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
