"""Regression tests for benchmark directory, resource, and case ID boundaries."""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pytest

from atst_tools.bench import batch_runner, batch_worker


def _options(tmp_path: Path, **overrides) -> batch_runner.BatchOptions:
    values = {"devices": ("0",), "output_dir": tmp_path / "out", "cpu_budget": 2,
              "telemetry": False, "case_telemetry": False}
    values.update(overrides)
    return batch_runner.BatchOptions(**values)


@pytest.mark.parametrize("alias", ["implicit", "explicit", "relative", "symlink"])
def test_resolved_workdirs_cannot_collide(tmp_path, alias):
    source = tmp_path / "source"
    source.mkdir()
    config_a = source / "a.yaml"
    config_b = source / "b.yaml"
    config_a.write_text("{}\n", encoding="utf-8")
    config_b.write_text("{}\n", encoding="utf-8")
    first = {"case_id": "a", "config": str(config_a)}
    second = {"case_id": "b", "config": str(config_b)}
    if alias == "explicit":
        second["workdir"] = str(source)
    elif alias == "relative":
        second["workdir"] = os.path.relpath(source, tmp_path / "out")
    elif alias == "symlink":
        link = tmp_path / "source-link"
        link.symlink_to(source, target_is_directory=True)
        second["workdir"] = str(link)
    with pytest.raises(ValueError, match="distinct directories"):
        batch_runner.load_cases({"cases": [first, second]}, _options(tmp_path))


def test_distinct_resolved_workdirs_remain_valid(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    rows = [{"case_id": "a", "config": str(source / "a.yaml")},
            {"case_id": "b", "config": str(source / "b.yaml"), "workdir": "b"}]
    assert len(batch_runner.load_cases({"cases": rows}, _options(tmp_path))) == 2


@pytest.mark.parametrize("case_id", ["../escape", "a/b", "a\\b", ".", "..", "bad\nline"])
def test_case_id_cannot_leave_its_report_directory(tmp_path, case_id):
    with pytest.raises(ValueError, match="case_id"):
        batch_runner.load_cases(
            {"cases": [{"case_id": case_id, "config": "input.yaml", "workdir": "safe"}]},
            _options(tmp_path),
        )
    with pytest.raises(ValueError, match="case_id"):
        batch_worker.Job.from_mapping(
            {"case_id": case_id, "config": "input.yaml", "workdir": str(tmp_path)}, 0
        )


def test_case_env_preserves_managed_facts_and_custom_values(tmp_path):
    overrides = {
        "CUDA_VISIBLE_DEVICES": "7", "OMP_NUM_THREADS": "99",
        "OPENBLAS_NUM_THREADS": "99", "ATST_THREADS_SOURCE": "explicit",
        "ATST_ATTEMPT": "99", "NUMBA_CACHE_DIR": "/tmp/shared",
        "ATST_RUNTIME_BOUND": "0", "ATST_INHERITED_DEVICES": "7",
        "ATST_EFFECTIVE_DEVICES": "7", "ATST_BINDING": "round_robin",
        "MY_FLAG": "kept",
    }
    case = batch_runner.CaseSpec.from_mapping(
        {"case_id": "one", "config": "input.yaml", "threads": 2, "env": overrides}
    )
    base = {"ATST_RUNTIME_BOUND": "1", "ATST_INHERITED_DEVICES": "0",
            "ATST_EFFECTIVE_DEVICES": "0", "ATST_BINDING": "inherit"}
    env = batch_runner.case_environment(case, ("0",), base=base, attempt=1, workdir=tmp_path)
    for key in batch_runner._launch.THREAD_ENV_KEYS:
        assert env[key] == "2"
    for key in batch_runner._launch.CACHE_ENV_KEYS:
        assert env[key] == str(tmp_path / ".atst_cache" / "attempt-1")
    assert env["CUDA_VISIBLE_DEVICES"] == "0"
    assert env["ATST_THREADS_SOURCE"] == "harness"
    assert env["ATST_ATTEMPT"] == "1"
    for key in ("ATST_RUNTIME_BOUND", "ATST_INHERITED_DEVICES", "ATST_EFFECTIVE_DEVICES", "ATST_BINDING"):
        assert env[key] == base[key]
    assert env["MY_FLAG"] == "kept"


@pytest.mark.parametrize("shared", [False, True])
def test_runtime_threads_and_runner_args_cannot_expand_batch_budget(tmp_path, shared):
    config = tmp_path / "input.yaml"
    config.write_text("runtime:\n  threads: 8\n", encoding="utf-8")
    options = _options(tmp_path, share_worker=shared)
    row = {"case_id": "one", "config": str(config), "workdir": "one", "threads": 2}
    with pytest.raises(ValueError, match="runtime.threads|thread budget"):
        batch_runner.run_manifest({"cases": [row]}, options)
    config.write_text("{}\n", encoding="utf-8")
    row["args"] = ["--threads", "8"]
    with pytest.raises(ValueError, match="runner option"):
        batch_runner.run_manifest({"cases": [row]}, options)


@pytest.mark.parametrize("shared", [False, True])
@pytest.mark.parametrize("backend", ["abacus", "dp", "legacy_abacus"])
@pytest.mark.parametrize("omp_value", ["8", "'8'"])
def test_calculator_explicit_omp_cannot_expand_batch_budget(tmp_path, shared, backend, omp_value):
    config = tmp_path / "input.yaml"
    prefix = "calculation:\n  type: relax\n  init_structure: missing.traj\n"
    if backend == "legacy_abacus":
        section = f"abacus:\n  omp: {omp_value}\n"
    else:
        required = "    model: missing.pb\n" if backend == "dp" else ""
        section = f"calculator:\n  name: {backend}\n  {backend}:\n{required}    omp: {omp_value}\n"
    config.write_text(prefix + section, encoding="utf-8")
    row = {"case_id": "one", "config": str(config), "workdir": "one", "threads": 2}
    field = "abacus" if backend == "legacy_abacus" else backend
    with pytest.raises(ValueError, match=f"calculator.{field}.omp"):
        batch_runner.run_manifest({"cases": [row]}, _options(tmp_path, share_worker=shared))
    config.write_text(prefix + section.replace(omp_value, "1"), encoding="utf-8")
    assert len(batch_runner.load_cases({"cases": [row]}, _options(tmp_path, share_worker=shared))) == 1


@pytest.mark.parametrize("option", ["--devices=7", "--dev=7", "--binding=round_robin", "--config=x.yaml",
                                     "--workdir=elsewhere", "--result-json=elsewhere.json"])
def test_runner_owned_args_cannot_change_resource_or_artifact_paths(tmp_path, option):
    row = {"case_id": "one", "config": "missing.yaml", "workdir": "one", "args": [option]}
    with pytest.raises(ValueError, match="runner option"):
        batch_runner.load_cases({"cases": [row]}, _options(tmp_path))


def test_runtime_auto_refused_when_affinity_exceeds_case_budget(tmp_path, monkeypatch):
    config = tmp_path / "input.yaml"
    config.write_text("runtime:\n  threads: auto\n", encoding="utf-8")
    monkeypatch.setattr(batch_runner._launch, "cpu_affinity_count", lambda: 8)
    row = {"case_id": "one", "config": str(config), "workdir": "one", "threads": 2}
    with pytest.raises(ValueError, match="runtime.threads"):
        batch_runner.load_cases({"cases": [row]}, _options(tmp_path))


@pytest.mark.parametrize("shared", [False, True])
def test_relative_config_identity_survives_worker_chdir(tmp_path, monkeypatch, shared):
    monkeypatch.chdir(tmp_path)
    config = tmp_path / "input.yaml"
    config.write_text("runtime:\n  threads: 1\n", encoding="utf-8")
    workdir = tmp_path / "out" / "case"
    workdir.mkdir(parents=True)
    (workdir / "input.yaml").write_text("runtime:\n  threads: 8\n", encoding="utf-8")
    script = tmp_path / "read_config.py"
    script.write_text(
        "import json,sys\nfrom pathlib import Path\n"
        "if len(sys.argv)==3:\n"
        " jobs=json.loads(Path(sys.argv[1]).read_text())\n"
        " for job in jobs:\n"
        "  source=Path(job['config'])\n"
        "  (Path(sys.argv[2])/'observed.txt').write_text(str(source.resolve())+'\\n'+source.read_text())\n"
        "  print(f\"[batch_worker] case={job['case_id']} status=success wall_s=0.001\",flush=True)\n"
        "else:\n"
        " source=Path(sys.argv[1])\n"
        " Path('observed.txt').write_text(str(source.resolve())+'\\n'+source.read_text())\n",
        encoding="utf-8",
    )
    options = _options(
        tmp_path, share_worker=shared,
        worker_factory=lambda case, cwd, devices: [sys.executable, str(script), case.config],
        shared_worker_factory=lambda jobs, out: [sys.executable, str(script), str(jobs), str(out)],
    )
    row = {"case_id": "case", "config": "input.yaml", "workdir": "case", "threads": 2}
    cases = batch_runner.load_cases({"cases": [row]}, options)
    summary = batch_runner.run_manifest({"cases": [row]}, options)
    assert summary["cases"][0]["status"] == "succeeded"
    observed = (tmp_path / "out" / "observed.txt" if shared else workdir / "observed.txt")
    assert observed.read_text(encoding="utf-8") == str(config.resolve()) + "\nruntime:\n  threads: 1\n"
    assert cases[0].config == str(config.resolve())
    assert batch_runner.worker_command(cases[0], workdir, ("0",))[4] == str(config.resolve())


@pytest.mark.parametrize("shared", [False, True])
def test_both_worker_modes_observe_managed_facts_without_runtime_yaml(tmp_path, shared):
    script = tmp_path / "facts_worker.py"
    script.write_text(
        "import json,os,sys\nfrom pathlib import Path\n"
        "keys=('CUDA_VISIBLE_DEVICES','OMP_NUM_THREADS','ATST_THREADS_SOURCE',"
        "'ATST_RUNTIME_BOUND','ATST_BINDING','NUMBA_CACHE_DIR','MY_FLAG')\n"
        "facts={key:os.environ.get(key) for key in keys}\n"
        "out=Path(sys.argv[1]) if len(sys.argv)>1 else Path.cwd()\n"
        "(out/'facts.json').write_text(json.dumps(facts))\n"
        "if len(sys.argv)>1:\n"
        " jobs=json.loads(Path(sys.argv[2]).read_text())\n"
        " for job in jobs:\n"
        "  print(f\"[batch_worker] case={job['case_id']} status=success wall_s=0.001\",flush=True)\n",
        encoding="utf-8",
    )
    env = {"CUDA_VISIBLE_DEVICES": "7", "OMP_NUM_THREADS": "99",
           "ATST_THREADS_SOURCE": "explicit", "ATST_RUNTIME_BOUND": "0",
           "ATST_BINDING": "round_robin", "NUMBA_CACHE_DIR": "/tmp/shared",
           "MY_FLAG": "kept"}
    options = _options(
        tmp_path, share_worker=shared,
        worker_factory=lambda case, workdir, devices: [sys.executable, str(script)],
        shared_worker_factory=lambda jobs, out: [sys.executable, str(script), str(out), str(jobs)],
    )
    row = {"case_id": "one", "config": "missing.yaml", "workdir": "one",
           "threads": 2, "env": env}
    summary = batch_runner.run_manifest({"cases": [row]}, options)
    assert summary["cases"][0]["status"] == "succeeded"
    facts_path = tmp_path / "out" / ("facts.json" if shared else "one/facts.json")
    facts = json.loads(facts_path.read_text(encoding="utf-8"))
    assert facts["CUDA_VISIBLE_DEVICES"] == "0"
    assert facts["OMP_NUM_THREADS"] == "2"
    assert facts["ATST_THREADS_SOURCE"] == "harness"
    assert facts["NUMBA_CACHE_DIR"].endswith("/.atst_cache/attempt-1")
    assert facts["MY_FLAG"] == "kept"
    if shared:
        assert facts["ATST_RUNTIME_BOUND"] == "1"
        assert facts["ATST_BINDING"] == "inherit"


@pytest.mark.parametrize("shared", [False, True])
def test_safe_space_id_runs_and_is_reported(tmp_path, shared):
    script = tmp_path / "worker.py"
    script.write_text(
        "import json,os,sys\nfrom pathlib import Path\n"
        "if len(sys.argv)==3:\n"
        " jobs=json.loads(Path(sys.argv[1]).read_text())\n"
        " for job in jobs:\n"
        "  p=Path(sys.argv[2])/job['case_id'];p.mkdir(parents=True,exist_ok=True)\n"
        "  (p/'atst_api_result.json').write_text('{}')\n"
        "  print(f\"[batch_worker] case={job['case_id']} status=success wall_s=0.001\",flush=True)\n"
        "else:\n"
        " (Path.cwd()/'atst_api_result.json').write_text('{}')\n",
        encoding="utf-8",
    )
    factory = lambda case, workdir, devices: [sys.executable, str(script)]
    shared_factory = lambda jobs, out: [sys.executable, str(script), str(jobs), str(out)]
    options = _options(tmp_path, share_worker=shared, worker_factory=factory,
                       shared_worker_factory=shared_factory)
    summary = batch_runner.run_manifest(
        {"cases": [{"case_id": "foo bar", "config": "input.yaml", "workdir": "safe"}]},
        options,
    )
    assert summary["cases"][0]["case_id"] == "foo bar"
    assert summary["cases"][0]["status"] == "succeeded"
    assert (tmp_path / "out" / "foo bar" / batch_runner.CASE_REPORT).is_file()
    assert batch_worker.parse_progress_line(
        batch_worker.format_progress_line("foo bar", "success", 0.001)
    ) == ("foo bar", "success", 0.001)
