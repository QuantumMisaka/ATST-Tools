"""Extract a normalized constant-potential record from an ATST result file."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any, Mapping


def extract_atst_point(
    results: Mapping[str, Any],
    point_index: int = 0,
) -> dict[str, Any]:
    """Convert one ``constant_potential_results.json`` point into a record.

    Args:
        results: Parsed ATST constant-potential result document.
        point_index: Zero-based index in ``results["points"]``.

    Returns:
        A JSON-compatible record accepted by
        :mod:`compare_constant_potential_results`.

    Raises:
        ValueError: If the document has no points or the index is out of range.
    """
    points = results.get("points")
    if not isinstance(points, list) or not points:
        raise ValueError("constant-potential results contain no points")
    try:
        point = points[point_index]
    except IndexError as exc:
        raise ValueError(
            f"point index {point_index} is outside 0..{len(points) - 1}"
        ) from exc
    if not isinstance(point, Mapping):
        raise ValueError(f"point {point_index} is not an object")

    record: dict[str, Any] = {}
    for source, target in (
        ("nelec", "nelec"),
        ("efermi", "efermi_eV"),
        ("vacuum_level", "vacuum_level_eV"),
        ("fermishift", "fermishift_eV"),
        ("mu_calc", "mu_calc_eV"),
        ("mu_target", "mu_target_eV"),
        ("residual_mu", "residual_mu_eV"),
        ("energy", "energy_eV"),
        ("free_energy", "free_energy_eV"),
        ("forces", "forces_eV_per_A"),
        ("cp_iterations", "iterations"),
        ("cp_converged", "converged"),
    ):
        if source in point:
            record[target] = point[source]

    compensation = point.get("compensation")
    if isinstance(compensation, Mapping):
        if "integrated_electrons" in compensation:
            record["density_integral_e"] = compensation["integrated_electrons"]
        if "compensation_derivative_ev" in compensation:
            record["compensation_derivative_eV"] = compensation[
                "compensation_derivative_ev"
            ]

    if point.get("capacitance_e_per_V") is not None:
        record["capacitance_e_per_V"] = point["capacitance_e_per_V"]
    else:
        history = point.get("cp_history")
        if isinstance(history, list):
            for item in reversed(history):
                if not isinstance(item, Mapping):
                    continue
                derivative = item.get("derivative_eV_per_e")
                try:
                    derivative_value = float(derivative)
                except (TypeError, ValueError):
                    continue
                if derivative_value > 0:
                    record["capacitance_e_per_V"] = 1.0 / derivative_value
                    break

    return record


def main(argv: list[str] | None = None) -> int:
    """Run the command-line extractor."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="ATST results JSON.")
    parser.add_argument("--output", required=True, help="Normalized record JSON.")
    parser.add_argument(
        "--point",
        type=int,
        default=0,
        help="Zero-based point index (default: 0).",
    )
    args = parser.parse_args(argv)

    results = json.loads(Path(args.input).read_text(encoding="utf-8"))
    try:
        record = extract_atst_point(results, args.point)
    except ValueError as exc:
        print(json.dumps({"ok": False, "error": str(exc)}, indent=2))
        return 2
    Path(args.output).write_text(
        json.dumps(record, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(record, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
