"""CPU-only confirmation of three ATST bench findings; writes under /tmp."""

from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

from atst_tools.bench import batch_runner, batch_worker
from atst_tools.runtime.cli_dispatch import plan_runner_launch
from atst_tools.runtime.launch import ensure_runtime_contract


def main() -> None:
    root = Path(tempfile.mkdtemp(prefix="atst-runtime-findings-", dir="/tmp"))
    print(f"fixture_root={root}")

    # 1. Distinct manifest spelling, identical actual work directory.
    output = root / "collision-out"
    config_dir = root / "same-workdir"
    options = batch_runner.BatchOptions(devices=("0",), output_dir=output)
    cases = batch_runner.load_cases(
        {
            "cases": [
                {"case_id": "default", "config": str(config_dir / "a.yaml")},
                {
                    "case_id": "explicit",
                    "config": str(root / "elsewhere" / "b.yaml"),
                    "workdir": str(config_dir),
                },
            ]
        },
        options,
    )
    paths = [batch_runner.case_workdir(case, output).resolve() for case in cases]
    print(f"workdir_collision.accepted={len(cases) == 2} actual_paths={paths}")
    assert paths[0] == paths[1]

    # 2. Real batch subprocess launch: custom env wins over allocated resources.
    probe = root / "env_probe.py"
    probe.write_text(
        "import json, os\n"
        "from pathlib import Path\n"
        "Path('observed-env.json').write_text(json.dumps({"
        "k: os.environ.get(k) for k in "
        "('CUDA_VISIBLE_DEVICES', 'OMP_NUM_THREADS', 'ATST_RUNTIME_BOUND')"
        "}))\n",
        encoding="utf-8",
    )
    (root / "env.yaml").write_text(
        "calculator:\n  name: abacus\ncalculation:\n  type: relax\n  init_structure: init.traj\n",
        encoding="utf-8",
    )
    env_output = root / "env-out"
    env_options = batch_runner.BatchOptions(
        devices=("0",),
        output_dir=env_output,
        cpu_budget=1,
        telemetry=False,
        case_telemetry=False,
        worker_factory=lambda case, workdir, devices: [sys.executable, str(probe)],
    )
    env_summary = batch_runner.run_manifest(
        {
            "cases": [
                {
                    "case_id": "env-case",
                    "config": str(root / "env.yaml"),
                    "workdir": "env-case",
                    "threads": 1,
                    "env": {
                        "CUDA_VISIBLE_DEVICES": "7",
                        "OMP_NUM_THREADS": "99",
                    },
                }
            ]
        },
        env_options,
    )
    observed = json.loads((env_output / "env-case" / "observed-env.json").read_text())
    print(
        "environment_override.status="
        f"{env_summary['cases'][0]['status']} allocated=0/1 observed={observed}"
    )
    assert observed["CUDA_VISIBLE_DEVICES"] == "7"
    assert observed["OMP_NUM_THREADS"] == "99"
    # A normal config with no runtime section takes the legacy runner route:
    # there is no later device/budget rebind or guard rejection.
    plan = plan_runner_launch(
        ["--config", str(root / "env.yaml"), "--workdir", str(env_output / "env-case")],
        environ={key: value for key, value in observed.items() if value is not None},
        workdir=env_output / "env-case",
    )
    ensure_runtime_contract(
        {}, environ={key: value for key, value in observed.items() if value is not None}
    )
    print(f"environment_override.runner_plan={plan!r} contract_guard=passed")
    assert plan is None

    # 3a. Case ID escapes output-dir report path.
    escaped = batch_runner.load_cases(
        {"cases": [{"case_id": "../escape", "config": "dummy.yaml", "workdir": "safe"}]},
        options,
    )[0]
    escaped_report = (output / escaped.case_id / batch_runner.CASE_REPORT).resolve()
    print(f"case_id_escape.accepted={escaped.case_id!r} report_path={escaped_report}")
    assert output.resolve() not in escaped_report.parents

    # 3b. Real shared-worker progress parser changes an accepted ID.
    shared_probe = root / "shared_probe.py"
    shared_probe.write_text(
        "import json, sys\n"
        "from pathlib import Path\n"
        "jobs = json.loads(Path(sys.argv[1]).read_text())\n"
        "for job in jobs:\n"
        "    p = Path(sys.argv[2]) / job['case_id']\n"
        "    p.mkdir(parents=True, exist_ok=True)\n"
        "    (p / 'atst_api_result.json').write_text('{}')\n"
        "    print(f\"[batch_worker] case={job['case_id']} status=success wall_s=0.001\", flush=True)\n",
        encoding="utf-8",
    )
    shared_output = root / "shared-out"
    shared_options = batch_runner.BatchOptions(
        devices=("0",),
        output_dir=shared_output,
        cpu_budget=1,
        telemetry=False,
        case_telemetry=False,
        share_worker=True,
        shared_worker_factory=lambda jobs, out: [sys.executable, str(shared_probe), str(jobs), str(out)],
    )
    shared_summary = batch_runner.run_manifest(
        {"cases": [{"case_id": "foo bar", "config": str(root / "shared.yaml"), "workdir": "shared"}]},
        shared_options,
    )
    progress = (shared_output / "worker.out").read_text().strip()
    parsed = batch_worker.parse_progress_line(progress)
    result_exists = (shared_output / "foo bar" / "atst_api_result.json").is_file()
    print(f"case_id_space.progress={progress!r} parsed={parsed!r}")
    print(
        "case_id_space.result_exists="
        f"{result_exists} worker_exit={shared_summary['shared_worker']['exit_code']} "
        f"report_status={shared_summary['cases'][0]['status']} "
        f"classification={shared_summary['cases'][0]['classification']}"
    )
    assert result_exists and parsed is not None and parsed[0] == "foo"
    assert shared_summary["cases"][0]["status"] == "failed"


if __name__ == "__main__":
    main()
