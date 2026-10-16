"""Workflow artifact manifest helpers."""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "atst-artifacts-v1"


def _jsonable(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if hasattr(value, "item"):
        return _jsonable(value.item())
    if isinstance(value, dict):
        return {key: _jsonable(val) for key, val in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def write_artifact_manifest(
    path: str | Path,
    *,
    workflow: str,
    artifacts: list[dict[str, Any]],
    stages: list[dict[str, Any]] | None = None,
    metadata: dict[str, Any] | None = None,
    plots: list[str] | tuple[str, ...] | None = None,
    backend_identity: dict[str, Any] | None = None,
    artifact_root: str | Path | None = None,
    path_base: str | Path | None = None,
) -> dict[str, Any]:
    """Write a stable machine-readable manifest for workflow outputs.

    Relative artifact paths are based at the writer's working directory, or
    ``artifact_root`` when provided. ``path_base`` preserves an existing
    manifest's explicit base during an enrichment rewrite. Both are expressed
    relative to the manifest directory on disk; root-level manifests omit the
    redundant ``.`` base for compatibility with existing v1 files.

    ``plots`` is an optional list of task-relative PNG paths produced for the
    workflow.  It is recorded under a dedicated top-level ``plots`` key only
    when non-empty, so manifests written without plots keep the exact original
    shape and older manifest consumers ignore the extension. ``backend_identity``
    is included only when explicit facts are provided by the caller.
    """
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "workflow": workflow,
        "metadata": metadata or {},
        "stages": stages or [],
        "artifacts": artifacts,
    }
    if plots:
        manifest["plots"] = [str(plot) for plot in plots]
    if backend_identity is not None:
        manifest["backend_identity"] = backend_identity
    output = Path(path)
    if artifact_root is not None and path_base is not None:
        raise ValueError("artifact_root and path_base cannot be used together")
    if path_base is None:
        root = Path.cwd() if artifact_root is None else Path(artifact_root).resolve()
        path_base = os.path.relpath(root, output.resolve().parent)
    if str(path_base) != ".":
        manifest["path_base"] = str(path_base)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(_jsonable(manifest), indent=2), encoding="utf-8")
    return manifest


def read_artifact_manifest(path: str | Path) -> dict[str, Any]:
    """Read an existing ATST artifact manifest without altering it."""
    return json.loads(Path(path).read_text(encoding="utf-8"))


def add_backend_identity_to_manifest(
    path: str | Path, identity: dict[str, Any]
) -> dict[str, Any]:
    """Add backend facts to an existing manifest without changing run facts.

    Args:
        path: Manifest already written by a workflow or the API finalizer.
        identity: Facts returned by ``calculator_backend_identity``.

    Returns:
        The enriched manifest as a detached dictionary.
    """
    output = Path(path)
    manifest = read_artifact_manifest(output)
    manifest["backend_identity"] = _jsonable(identity)
    output.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest
