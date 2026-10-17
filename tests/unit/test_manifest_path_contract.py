"""Artifact path handoff between workflow writers, summaries, and the API."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from helpers import FakeWorld

from atst_tools.api import RunOptions, run_workflow
from atst_tools.api import services
from atst_tools.utils.artifacts import (
    add_backend_identity_to_manifest,
    read_artifact_manifest,
    write_artifact_manifest,
)
from atst_tools.utils.summary import summarize_artifact_manifest


def test_nested_manifest_resolves_workdir_paths_from_another_cwd(monkeypatch, tmp_path):
    workdir = tmp_path / "run"
    workdir.mkdir()
    trajectory = workdir / "neb.traj"
    trajectory.write_text("frames", encoding="utf-8")
    absolute = workdir / "absolute.traj"
    absolute.write_text("frames", encoding="utf-8")
    monkeypatch.chdir(workdir)

    written = write_artifact_manifest(
        "reports/nested/artifacts.json",
        workflow="neb",
        artifacts=[
            {"role": "trajectory", "path": "neb.traj"},
            {"role": "absolute", "path": str(absolute)},
            {"role": "missing", "path": "missing.traj"},
        ],
    )
    manifest_path = workdir / "reports/nested/artifacts.json"
    monkeypatch.chdir(tmp_path)
    summary = summarize_artifact_manifest(manifest_path)

    assert written["path_base"] == "../.."
    assert [
        artifact["path"]
        for artifact in read_artifact_manifest(manifest_path)["artifacts"]
    ] == [
        "neb.traj", str(absolute), "missing.traj"
    ]
    assert [item["exists"] for item in summary["artifacts"]] == [True, True, False]
    assert Path(summary["artifacts"][0]["resolved_path"]).resolve() == trajectory
    assert Path(summary["artifacts"][1]["resolved_path"]) == absolute
    assert (
        Path(summary["artifacts"][2]["resolved_path"]).resolve()
        == workdir / "missing.traj"
    )
    assert summary["status"]["n_missing"] == 1


def test_legacy_root_manifest_keeps_manifest_relative_paths(monkeypatch, tmp_path):
    workdir = tmp_path / "run"
    workdir.mkdir()
    (workdir / "neb.traj").write_text("frames", encoding="utf-8")
    (workdir / "atst_artifacts.json").write_text(
        json.dumps(
            {
                "schema_version": "atst-artifacts-v1",
                "workflow": "neb",
                "artifacts": [{"path": "neb.traj"}],
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.chdir(tmp_path)

    summary = summarize_artifact_manifest(workdir / "atst_artifacts.json")

    assert summary["status"]["complete"] is True
    assert Path(summary["artifacts"][0]["resolved_path"]) == workdir / "neb.traj"

    monkeypatch.chdir(workdir)
    root_manifest = write_artifact_manifest(
        "new_artifacts.json",
        workflow="neb",
        artifacts=[{"path": "neb.traj"}],
    )
    assert "path_base" not in root_manifest


@pytest.mark.parametrize("external_manifest", [False, True])
def test_api_rewrite_preserves_nested_manifest_path_base(
    monkeypatch, tmp_path, external_manifest
):
    workdir = tmp_path / "run"
    workdir.mkdir()
    monkeypatch.chdir(workdir)
    (workdir / "neb.traj").write_text("frames", encoding="utf-8")
    manifest_name = (
        str(tmp_path / "external/artifacts.json")
        if external_manifest
        else "reports/artifacts.json"
    )

    def complete_neb(config, options):
        write_artifact_manifest(
            config["calculation"]["artifact_manifest"],
            workflow="neb",
            artifacts=[{"role": "trajectory", "path": "neb.traj"}],
            metadata={"manifest_source": "runner"},
        )
        return None

    monkeypatch.setattr(services, "_dispatch_normalized", complete_neb)
    monkeypatch.setattr(services, "_workflow_plots", lambda *args: ("neb_profile.png",))
    config = {
        "calculation": {
            "type": "neb",
            "init_chain": "chain.traj",
            "artifact_manifest": manifest_name,
        },
        "calculator": {"name": "abacus", "abacus": {"parameters": {}}},
    }

    result = run_workflow(config, RunOptions(world=FakeWorld(), plots=True))
    document = result.to_document(workdir)
    manifest_path = Path(document["artifact_manifest"])
    first = read_artifact_manifest(manifest_path)
    add_backend_identity_to_manifest(manifest_path, first["backend_identity"])
    assert read_artifact_manifest(manifest_path) == first
    monkeypatch.chdir(tmp_path.parent)
    summary = summarize_artifact_manifest(manifest_path)

    assert first["path_base"] == ("../run" if external_manifest else "..")
    assert first["artifacts"] == [{"role": "trajectory", "path": "neb.traj"}]
    assert first["plots"] == ["neb_profile.png"]
    assert first["metadata"]["manifest_source"] == "runner"
    assert document["artifacts"] == first["artifacts"]
    assert summary["status"]["complete"] is True
    assert (
        Path(summary["artifacts"][0]["resolved_path"]).resolve()
        == workdir / "neb.traj"
    )


def test_api_synthesized_manifest_uses_workdir_base(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "neb.traj").write_text("frames", encoding="utf-8")
    monkeypatch.setattr(services, "_dispatch_normalized", lambda *args: None)
    config = {
        "calculation": {
            "type": "neb",
            "init_chain": "chain.traj",
            "artifact_manifest": "reports/artifacts.json",
        },
        "calculator": {"name": "abacus", "abacus": {"parameters": {}}},
    }

    result = run_workflow(config, RunOptions(world=FakeWorld()))
    manifest_path = Path(result.to_document(tmp_path)["artifact_manifest"])
    summary = summarize_artifact_manifest(manifest_path)

    assert read_artifact_manifest(manifest_path)["path_base"] == ".."
    assert result.artifacts == ({"role": "trajectory", "path": "neb.traj"},)
    assert summary["status"]["complete"] is True
