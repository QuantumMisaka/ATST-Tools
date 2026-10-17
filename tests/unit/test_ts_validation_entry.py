"""Public entry points for explicit stored-data TS screening."""

from __future__ import annotations

import json

import pytest


def _vibration_facts():
    geometry = {
        "numbers": [1],
        "positions": [[0.0, 0.0, 0.0]],
        "cell": [[0.0, 0.0, 0.0]] * 3,
        "pbc": [False, False, False],
    }
    return {"imaginary_frequencies": [120.0], "n_atoms": 1, "indices": [0], "geometry": geometry}


def test_public_api_exposes_explicit_validation_without_running_calculator(tmp_path):
    from atst_tools.api import validate_ts

    vibration = tmp_path / "vibration.json"
    vibration.write_text(json.dumps(_vibration_facts()), encoding="utf-8")

    summary = validate_ts(vibration)

    assert summary["source"] == str(vibration)
    assert summary["status"] == "pass"
    assert summary["stationarity"]["status"] == "unknown"
    assert summary["requested_checks_status"] == "unknown"
    assert sorted(path.name for path in tmp_path.iterdir()) == ["vibration.json"]


def test_cli_validation_reads_existing_results_and_writes_only_requested_output(tmp_path, capsys):
    from atst_tools.scripts.cli import main

    vibration = tmp_path / "vibration.json"
    structure = tmp_path / "structure.json"
    output = tmp_path / "validation.json"
    results = _vibration_facts()
    vibration.write_text(json.dumps(results), encoding="utf-8")
    structure.write_text(
        json.dumps(
            {
                "n_atoms": 1,
                "geometry": results["geometry"],
                "fixed_indices": [],
                "forces": [[0.03, 0.0, 0.0]],
            }
        ),
        encoding="utf-8",
    )

    assert main([
        "vibration", "validate", str(vibration), "--structure", str(structure),
        "--fmax-threshold", "0.02", "--output", str(output),
    ]) is None

    summary = json.loads(capsys.readouterr().out)
    assert summary == json.loads(output.read_text(encoding="utf-8"))
    assert summary["stationarity"]["status"] == "fail"
    assert summary["requested_checks_status"] == "fail"
    assert json.loads(vibration.read_text(encoding="utf-8")) == results


def test_cli_validation_without_structure_reports_unknown_and_creates_no_output(tmp_path, capsys):
    from atst_tools.scripts.cli import main

    vibration = tmp_path / "vibration.json"
    vibration.write_text(json.dumps(_vibration_facts()), encoding="utf-8")

    main(["vibration", "validate", str(vibration)])

    assert json.loads(capsys.readouterr().out)["requested_checks_status"] == "unknown"
    assert sorted(path.name for path in tmp_path.iterdir()) == ["vibration.json"]


def test_validation_cannot_overwrite_either_input(tmp_path):
    from atst_tools.api import validate_ts

    vibration = tmp_path / "vibration.json"
    structure = tmp_path / "structure.json"
    results = _vibration_facts()
    vibration.write_text(json.dumps(results), encoding="utf-8")
    structure.write_text(json.dumps({"n_atoms": 1, "geometry": results["geometry"]}), encoding="utf-8")
    original_vibration = vibration.read_bytes()
    original_structure = structure.read_text(encoding="utf-8")
    alias = tmp_path / "structure_alias.json"
    alias.symlink_to(structure)

    with pytest.raises(ValueError, match="must differ from.*input"):
        validate_ts(vibration, output=vibration)
    with pytest.raises(ValueError, match="must differ from.*input"):
        validate_ts(vibration, structure=structure, output=structure)
    with pytest.raises(ValueError, match="must differ from.*input"):
        validate_ts(vibration, structure=structure, output=alias)

    assert vibration.read_bytes() == original_vibration
    assert structure.read_text(encoding="utf-8") == original_structure
