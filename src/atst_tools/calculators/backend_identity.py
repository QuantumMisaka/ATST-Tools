"""Lightweight, durable facts about the calculator selected for a run."""

from __future__ import annotations

import hashlib
import re
from importlib import metadata
from pathlib import Path
from typing import Any, Mapping

from atst_tools.calculators.abacuslite_backend import Abacus, BACKEND_SOURCE
from atst_tools.calculators.dp import DP_BACKEND_BY_KIND, dp_model_kind, dp_section


_VENDORED_ROOT = Path(__file__).resolve().parents[1] / "external" / "ASE_interface"
_SNAPSHOT_DOCUMENT = "atst_tools/external/ASE_interface/ABACUSLITE_SNAPSHOT.md"
_PATCHES_DOCUMENT = "atst_tools/external/ASE_interface/PATCHES.md"


def _package_version(distribution: str) -> str:
    """Return installed distribution metadata without importing its runtime."""
    try:
        return metadata.version(distribution)
    except metadata.PackageNotFoundError:
        return "unknown"


def _calculator_type(calculator: Any) -> str:
    """Return the fully qualified type of an existing calculator object."""
    cls = type(calculator)
    return f"{cls.__module__}.{cls.__qualname__}"


def _vendored_provenance() -> dict[str, str] | str:
    """Read the packaged snapshot marker and small patch ledger."""
    snapshot_path = _VENDORED_ROOT / "ABACUSLITE_SNAPSHOT.md"
    patches_path = _VENDORED_ROOT / "PATCHES.md"
    try:
        snapshot_text = snapshot_path.read_text(encoding="utf-8")
        patches = patches_path.read_bytes()
    except OSError:
        return "unknown"
    match = re.search(r"当前基线：上游\s+`([0-9a-f]{40})`", snapshot_text)
    if match is None:
        return "unknown"
    return {
        "upstream_sha": match.group(1),
        "snapshot_document": _SNAPSHOT_DOCUMENT,
        "patches_document": _PATCHES_DOCUMENT,
        "patches_sha256": hashlib.sha256(patches).hexdigest(),
    }


def _model_facts(model: Any, head: Any) -> dict[str, Any] | str:
    """Describe an already selected DP model path without opening the model."""
    if not model:
        return "unknown"
    path = Path(str(model)).expanduser().resolve()
    kind = dp_model_kind(path)
    try:
        size: int | str = path.stat().st_size
    except OSError:
        size = "unknown"
    return {
        "path": str(path),
        "kind": kind,
        "backend": DP_BACKEND_BY_KIND.get(kind) or "unknown",
        "head": str(head) if head is not None else "unknown",
        "size_bytes": size,
    }


def calculator_backend_identity(
    config: Mapping[str, Any] | None = None,
    *,
    calculator: Any | None = None,
    supplied: bool = False,
) -> dict[str, Any]:
    """Record available backend facts without constructing or evaluating it.

    Args:
        config: Normalized configuration for a factory-owned calculator.
        calculator: Existing calculator instance, if one is available.
        supplied: Whether the caller owns the calculator rather than the factory.

    Returns:
        JSON-compatible identity with an explicit fact basis. A configuration
        identity records the resolved backend selection; only an instance
        identity claims that the calculator object was observed. Unavailable
        versions, models and provenance are written as ``"unknown"``.
    """
    if supplied:
        return {
            "calculator": "supplied",
            "identity_basis": "calculator_instance" if calculator is not None else "configuration",
            "source": "caller",
            "calculator_type": _calculator_type(calculator) if calculator is not None else "unknown",
            "backend_version": "unknown",
            "snapshot": "unknown",
            "model": _model_facts(
                getattr(calculator, "atst_model", None),
                getattr(calculator, "atst_model_head", None),
            ) if calculator is not None else "unknown",
        }

    configured = config or {}
    section = configured.get("calculator", {})
    section = section if isinstance(section, Mapping) else {}
    name = str(section.get("name", "abacus")).lower()
    basis = "calculator_instance" if calculator is not None else "configuration"
    type_name = _calculator_type(calculator) if calculator is not None else "unknown"

    if name == "abacus":
        source = BACKEND_SOURCE if calculator is None or isinstance(calculator, Abacus) else "unknown"
        return {
            "calculator": "abacus",
            "identity_basis": basis,
            "source": source,
            "calculator_type": type_name,
            "backend_version": _package_version("abacuslite") if source == "external" else "unknown",
            "snapshot": _vendored_provenance() if source == "vendored" else "unknown",
            "model": "unknown",
        }

    if name in {"dp", "deepmd"}:
        # The factory attaches these attributes to its completed calculator.
        # A caller-supplied object without them cannot establish DP provenance.
        factory_instance = calculator is not None and hasattr(calculator, "atst_model")
        source = "deepmd" if calculator is None or factory_instance else "unknown"
        dp_config = dp_section(dict(configured))
        if factory_instance:
            model = getattr(calculator, "atst_model", None)
            head = getattr(calculator, "atst_model_head", None)
        elif calculator is None:
            model = dp_config.get("model")
            head = dp_config.get("head")
        else:
            model = head = None
        return {
            "calculator": "dp",
            "identity_basis": basis,
            "source": source,
            "calculator_type": type_name,
            "backend_version": _package_version("deepmd-kit") if source == "deepmd" else "unknown",
            "snapshot": "unknown",
            "model": _model_facts(model, head),
        }

    return {
        "calculator": name,
        "identity_basis": basis,
        "source": "unknown",
        "calculator_type": type_name,
        "backend_version": "unknown",
        "snapshot": "unknown",
        "model": "unknown",
    }
