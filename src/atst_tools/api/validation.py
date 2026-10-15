"""Explicit read-only transition-state screening entry point."""

from __future__ import annotations

import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from atst_tools.utils.ts_validation import build_ts_validation_summary, frozen_structure_fact


def _read_mapping(source: Mapping[str, Any] | str | Path, name: str) -> tuple[dict[str, Any], str | None]:
    """Load a JSON fact document or copy a supplied mapping."""
    if isinstance(source, Mapping):
        return dict(source), None
    if not isinstance(source, (str, Path)):
        raise TypeError(f"{name} must be a mapping or path.")
    path = Path(source)
    contents = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(contents, dict):
        raise TypeError(f"{name} JSON must contain an object.")
    return contents, str(path)


def validate_ts(
    vibration_results: Mapping[str, Any] | str | Path,
    *,
    structure: Mapping[str, Any] | str | Path | None = None,
    fmax_threshold: float = 0.05,
    imaginary_cutoff_cm1: float = 50.0,
    output: str | Path | None = None,
) -> dict[str, Any]:
    """Screen existing vibration output against optional frozen structure facts.

    This operation reads stored data only. A structure JSON may contain paired
    geometry, atom count, fixed indices and forces. An ASE structure path may
    contain frozen ``SinglePointCalculator`` forces; a live calculator is never
    queried. The result is scoped screening, not a proof of IRC connectivity.

    Args:
        vibration_results: Vibration JSON path or already loaded mapping.
        structure: Paired frozen evidence JSON mapping/path or ASE structure
            path. Its geometry must match that in the vibration results before
            force evidence is accepted.
        fmax_threshold: Stationarity threshold in eV/Angstrom.
        imaginary_cutoff_cm1: Significant imaginary frequency cutoff in cm^-1.
        output: Optional path for the validation summary JSON.

    Returns:
        JSON-ready validation summary.

    Raises:
        ValueError: Supplied facts are malformed, nonfinite or inconsistent.
        TypeError: A supplied document has the wrong type.
    """
    results, source = _read_mapping(vibration_results, "vibration_results")
    fact = None
    if structure is not None:
        if isinstance(structure, (str, Path)) and Path(structure).suffix.lower() != ".json":
            from ase.io import read

            fact = frozen_structure_fact(read(structure, index=-1))
        else:
            fact, _ = _read_mapping(structure, "structure")
    summary = build_ts_validation_summary(
        results,
        fmax_threshold=fmax_threshold,
        imaginary_cutoff_cm1=imaginary_cutoff_cm1,
        source=source,
        structure=fact,
    )
    if output is not None:
        output_path = Path(output)
        if source is not None and output_path.resolve() == Path(source).resolve():
            raise ValueError("Validation output must differ from vibration input.")
        if isinstance(structure, (str, Path)) and output_path.resolve() == Path(structure).resolve():
            raise ValueError("Validation output must differ from structure input.")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary
