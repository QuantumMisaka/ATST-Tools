"""Probe existing platform MPI entry contracts without executing ABACUS.

Run with abacus-env and pass the parent ABACUS toolbox root. External runtime
refresh, scientific processes, and plotting are stubbed; mapping validation,
configuration construction, resource hooks, and the runner boundary are real.
This is an audit reproducer for the recorded baseline, not a regression suite.
"""

from __future__ import annotations

import contextlib
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import pytest
import yaml


def main() -> None:
    """Print bounded evidence for prepare/run allocation inconsistencies."""
    root = Path(sys.argv[1]).resolve()
    sys.path[:0] = [str(root), str(root / "tests")]
    from adam_test_utils import install_fake_runner, install_tool_stub

    files = ["abacus_transition.py", "toolkits/transition_atst.py", "utils/env_setting.py"]
    result = {
        "scope": "offline stubs; no MPI, GPU, ABACUS, scheduler, or platform execution",
        "parent_head": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=root, text=True
        ).strip(),
        "atst_head": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=root / "deps/atst-tools", text=True
        ).strip(),
        "source_sha256": {
            name: hashlib.sha256((root / name).read_bytes()).hexdigest() for name in files
        },
        "python": sys.executable,
    }
    with pytest.MonkeyPatch.context() as patch, tempfile.TemporaryDirectory() as folder:
        install_tool_stub(patch)
        import abacus_transition as tools
        import toolkits.transition_atst as transition
        import utils.env_setting as runtime

        prepare = tools.abacus_worker_transition_prepare_neb()
        run = tools.abacus_worker_transition_run()
        result["allocation"] = {
            "prepare_gpu4_intent": prepare.resAlloc({"args": {"n_gpu": 4}, "tip": ""}),
            "run_explicit_gpu4": run.resAlloc({"args": {"n_gpu": 4}, "tip": ""}),
            "run_omitted_gpu": run.resAlloc({"args": {}, "tip": ""}),
        }
        # Preserve the real shape/mapping gate. Only SIF refresh is disabled.
        patch.setattr(runtime, "refresh_sai_runtime_from_sif", lambda **kwargs: None)
        patch.setattr(runtime, "abacus_runtime_ready", lambda: True)
        for key in list(runtime.os.environ):
            if key.startswith(("SLURM_", "SAI_", "ABACUS_", "OMPI_", "PMI", "ATST_")):
                patch.delenv(key, raising=False)
        for key, value in {
            "SAI_NNODES": "1", "SAI_GPUS_ON_NODE": "1", "SAI_NTASKS": "1",
            "SAI_MAPPING_STATUS": "loaded", "RANKS_PER_GPU": "1",
            "MAP_OPT": "ppr:1:node:pe=4",
        }.items():
            patch.setenv(key, value)
        args = dict(
            structure_paths=[], n_cpu=16, n_gpu=4, mpi_mode=1, kpts=[1, 1, 1],
            abacus_parameters={"ks_solver": "cusolver"},
            pseudopotentials={"H": "H.upf"}, basissets={"H": "H.orb"},
            pseudo_dir=".", orbital_dir=".",
        )
        with contextlib.redirect_stdout(io.StringIO()):
            try:
                transition._build_abacus_calculator(**args)
            except runtime.SaiRuntimeContractError as exc:
                result["prepare_gpu4_on_gpu1"] = {"rejected": True, "error": str(exc)}
            else:
                result["prepare_gpu4_on_gpu1"] = {"rejected": False}
            for key, value in {
                "SAI_GPUS_ON_NODE": "4", "SAI_NTASKS": "4",
                "SAI_MAPPING_STATUS": "loaded", "MAP_OPT": "ppr:4:node:pe=4",
            }.items():
                patch.setenv(key, value)
            calculator = transition._build_abacus_calculator(**args)
        result["calculator_on_gpu4"] = {
            key: calculator["abacus"][key] for key in ("command", "mpi", "omp")
        }
        try:
            transition._validate_image_parallel_calculator({"calculator": calculator})
        except ValueError as exc:
            result["nested_mpi"] = {"rejected": True, "error": str(exc)}
        else:
            result["nested_mpi"] = {"rejected": False}

        # Execute the real run adapter up to its mocked subprocess boundary.
        # The saved command requests four MPI ranks; the caller omits n_gpu.
        work = Path(folder)
        config = work / "atst_neb.yaml"
        config.write_text(yaml.safe_dump({
            "calculation": {"type": "neb", "parallel": False, "init_chain": "chain.traj"},
            "calculator": calculator,
        }), encoding="utf-8")
        calls = install_fake_runner(patch, transition, record="cmd")
        patch.setattr(transition, "render_transition_plots", lambda **kwargs: [])
        patch.setattr(transition, "_is_atst_check_input_runtime_unavailable", lambda: False)
        with contextlib.redirect_stdout(io.StringIO()):
            outcome = transition.run_transition_config(config_path=str(config))
        result["run_gpu4_command_without_gpu_arg"] = {
            "adapter_status_with_stub_runner": outcome.status,
            "formal_runner_reached": any("--dry-run" not in call for call in calls),
            "saved_command": outcome.normalized_config["calculator"]["abacus"]["command"],
            "allocated_gpus": result["allocation"]["run_omitted_gpu"]["GPU"],
            "outer_image_parallel": outcome.outer_command_evidence["image_parallel"],
        }
        mapping_requests = []
        patch.setattr(transition, "ensure_sai_runtime_mapping", lambda **kwargs: mapping_requests.append(kwargs))
        patch.setenv("MAP_OPT", "ppr:2:package:pe=4")
        with contextlib.redirect_stdout(io.StringIO()):
            explicit = transition.run_transition_config(config_path=str(config), n_gpu=4)
        result["run_gpu4_with_changed_mapping"] = {
            "adapter_status_with_stub_runner": explicit.status,
            "mapping_refresh_calls": mapping_requests,
            "current_map_opt": "ppr:2:package:pe=4",
            "saved_command": explicit.normalized_config["calculator"]["abacus"]["command"],
        }
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
