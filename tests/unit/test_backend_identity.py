"""Backend identity records facts without launching or reloading calculators."""

from __future__ import annotations

import json
from pathlib import Path

from atst_tools.calculators.backend_identity import calculator_backend_identity
from atst_tools.utils.artifacts import (
    add_backend_identity_to_manifest,
    write_artifact_manifest,
)


def test_configured_abacus_records_external_source_and_version():
    identity = calculator_backend_identity({"calculator": {"name": "abacus"}})

    assert identity["calculator"] == "abacus"
    assert identity["identity_basis"] == "configuration"
    assert identity["source"] == "external"
    assert identity["backend_version"] != ""
    assert identity["snapshot"] == "unknown"


def test_configured_dp_uses_model_and_head_facts_without_loading_model(tmp_path, monkeypatch):
    model = tmp_path / "model.pth"
    model.write_bytes(b"placeholder")
    monkeypatch.setattr(
        "atst_tools.calculators.dp.dp_model_identity",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("model reload")),
    )

    identity = calculator_backend_identity(
        {"calculator": {"name": "dp", "dp": {"model": str(model), "head": "branch"}}}
    )

    assert identity["calculator"] == "dp"
    assert identity["identity_basis"] == "configuration"
    assert identity["source"] == "deepmd"
    assert identity["model"] == {
        "path": str(model),
        "kind": "frozen_archive",
        "backend": "pt",
        "head": "branch",
        "size_bytes": len(b"placeholder"),
    }


def test_supplied_calculator_records_type_without_evaluation():
    class CallerCalculator:
        def get_potential_energy(self):
            raise AssertionError("unexpected evaluation")

    calculator = CallerCalculator()
    identity = calculator_backend_identity(calculator=calculator, supplied=True)

    assert identity["calculator"] == "supplied"
    assert identity["identity_basis"] == "calculator_instance"
    assert identity["source"] == "caller"
    assert identity["calculator_type"].endswith("CallerCalculator")
    assert identity["backend_version"] == "unknown"
    assert identity["model"] == "unknown"


def test_dp_instance_without_factory_facts_keeps_source_and_model_unknown(tmp_path):
    model = tmp_path / "configured.pb"
    model.write_bytes(b"placeholder")
    identity = calculator_backend_identity(
        {"calculator": {"name": "dp", "dp": {"model": str(model)}}},
        calculator=object(),
    )

    assert identity["identity_basis"] == "calculator_instance"
    assert identity["source"] == "unknown"
    assert identity["model"] == "unknown"


def test_sequential_configurations_do_not_reuse_previous_model(tmp_path):
    first = tmp_path / "first.pt"
    second = tmp_path / "second.pb"
    first.write_bytes(b"first")
    second.write_bytes(b"second")

    first_identity = calculator_backend_identity(
        {"calculator": {"name": "dp", "dp": {"model": str(first), "head": "A"}}}
    )
    second_identity = calculator_backend_identity(
        {"calculator": {"name": "dp", "dp": {"model": str(second)}}}
    )

    assert first_identity["model"]["path"] == str(first)
    assert first_identity["model"]["head"] == "A"
    assert second_identity["model"]["path"] == str(second)
    assert second_identity["model"]["head"] == "unknown"


def test_manifest_enrichment_preserves_workflow_facts(tmp_path):
    path = tmp_path / "artifacts.json"
    original = write_artifact_manifest(
        path,
        workflow="ccqn",
        artifacts=[{"role": "trajectory", "path": "ccqn.traj"}],
        stages=[{"name": "ccqn", "status": "complete", "converged": None}],
        metadata={"existing": True},
    )
    identity = {"calculator": "supplied", "source": "caller"}

    enriched = add_backend_identity_to_manifest(path, identity)

    assert enriched["backend_identity"] == identity
    assert {key: value for key, value in enriched.items() if key != "backend_identity"} == original
    assert json.loads(path.read_text(encoding="utf-8")) == enriched
