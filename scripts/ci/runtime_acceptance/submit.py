#!/usr/bin/env python3
"""Package installed wheels and stream one bounded SAI runtime acceptance job."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import tarfile
from datetime import datetime, timezone

HERE = Path(__file__).resolve().parent
REQUIRED_ASSETS = (
    "h2.stru", "pt.stru", "h.upf", "h.orb", "pt.upf", "pt.orb",
    "h2_chain.traj", "h2_mode.npy",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def required_file(value: str | Path) -> Path:
    path = Path(value).expanduser().resolve(strict=True)
    if not path.is_file() or not path.stat().st_size:
        raise ValueError(f"not a nonempty file: {path}")
    return path


def full_sha(value: str, name: str) -> str:
    if len(value) != 40 or any(char not in "0123456789abcdef" for char in value):
        raise ValueError(f"{name} must be a lowercase full Git SHA")
    return value


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atst-wheel", required=True)
    parser.add_argument("--atst-source-sha", required=True)
    parser.add_argument("--plugin-wheel", required=True)
    parser.add_argument("--plugin-source-sha", required=True)
    parser.add_argument("--engine", required=True)
    parser.add_argument("--engine-sha256", required=True)
    parser.add_argument("--engine-source-sha", required=True)
    parser.add_argument("--sella-wheel", required=True)
    parser.add_argument("--assets", required=True, help="directory with the eight named input assets")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--execute", action="store_true", help="submit the prepared kit to SAI")
    args = parser.parse_args()
    output = Path(args.output_dir).expanduser().resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"output directory must be new or empty: {output}")
    output.mkdir(parents=True, exist_ok=True)
    atst = required_file(args.atst_wheel)
    plugin = required_file(args.plugin_wheel)
    engine = required_file(args.engine)
    sella = required_file(args.sella_wheel)
    assets = Path(args.assets).expanduser().resolve(strict=True)
    if not assets.is_dir():
        raise ValueError(f"not an assets directory: {assets}")
    asset_files = {name: required_file(assets / name) for name in REQUIRED_ASSETS}
    if not atst.name.startswith("atst_tools-") or atst.suffix != ".whl":
        raise ValueError(f"unexpected ATST wheel name: {atst.name}")
    if not plugin.name.startswith("abacuslite-") or plugin.suffix != ".whl":
        raise ValueError(f"unexpected plugin wheel name: {plugin.name}")
    if sha256(engine) != args.engine_sha256:
        raise ValueError("engine bytes differ from expected SHA256")
    atst_sha = full_sha(args.atst_source_sha, "ATST source")
    plugin_sha = full_sha(args.plugin_source_sha, "plugin source")
    engine_sha = full_sha(args.engine_source_sha, "engine source")
    wheels = {"atst": atst, "plugin": plugin, "sella": sella}
    identity = {
        "schema": "atst-runtime-acceptance-source-v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "atst": {"source_sha": atst_sha, "wheel": str(atst), "sha256": sha256(atst)},
        "abacuslite": {"source_sha": plugin_sha, "wheel": str(plugin), "sha256": sha256(plugin)},
        "engine": {"source_sha": engine_sha, "binary": str(engine), "sha256": sha256(engine)},
        "sella": {"wheel": str(sella), "sha256": sha256(sella)},
        "assets": {name: {"source": str(path), "sha256": sha256(path)}
                   for name, path in asset_files.items()},
        "harness": {name: sha256(HERE / name) for name in ("acceptance.py", "run.sh", "submit.py")},
        "scope": "real-engine interface and output acceptance; no scientific convergence claim",
    }
    identity_file = output / "source-identity.json"
    identity_file.write_text(json.dumps(identity, indent=2, sort_keys=True) + "\n")
    pins = {
        "EXPECTED_ATST_SHA256": identity["atst"]["sha256"],
        "EXPECTED_PLUGIN_SHA256": identity["abacuslite"]["sha256"],
        "ATST_SOURCE_SHA": atst_sha,
        "PLUGIN_SOURCE_SHA": plugin_sha,
        "ABACUS_SOURCE_SHA": engine_sha,
        "ATST_WHEEL_NAME": atst.name,
        "PLUGIN_WHEEL_NAME": plugin.name,
        "SELLA_WHEEL_NAME": sella.name,
    }
    pin_file = output / "pins.sh"
    pin_file.write_text("".join(f"export {key}={shlex.quote(value)}\n" for key, value in pins.items()))
    kit = output / "runtime-payload.tar"
    with tarfile.open(kit, "w") as archive:
        for path, name in ((HERE / "run.sh", "run.sh"),
                           (HERE / "acceptance.py", "harness/acceptance.py"),
                           (identity_file, "source-identity.json"),
                           (pin_file, "pins.sh"), (engine, "abacus_basic_gpu")):
            archive.add(path, arcname=name)
        for name, path in asset_files.items():
            archive.add(path, arcname="harness/assets/" + name)
        for path in wheels.values():
            archive.add(path, arcname="wheels/" + path.name)
    kit_meta = {"kit": str(kit), "sha256": sha256(kit), "bytes": kit.stat().st_size,
                "identity": str(identity_file), "execute": args.execute}
    (output / "kit.json").write_text(json.dumps(kit_meta, indent=2, sort_keys=True) + "\n")
    print(json.dumps(kit_meta, sort_keys=True), flush=True)
    if not args.execute:
        return
    entry = (
        'set -e; task_root=$(mktemp -d /tmp/atst-runtime.XXXXXX); '
        'tar -xf - -C "$task_root"; exec bash "$task_root/run.sh" "$task_root"'
    )
    remote = "cd /tmp && " + shlex.join([
        "srun", "--account=pku-jianghong", "--partition=4V100", "--qos=improper-gpu",
        "--nodes=1", "--ntasks=1", "--gpus-per-node=2", "--time=03:00:00",
        "--job-name=atst-runtime-acceptance", "bash", "-c", entry,
    ])
    result_path = output / "runtime-result.tar.gz"
    log_path = output / "runtime-sai.log"
    with kit.open("rb") as source, result_path.open("wb") as result, log_path.open("wb") as log:
        status = subprocess.run([
            "ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=15",
            "-o", "ServerAliveInterval=30", "SAI", remote,
        ], stdin=source, stdout=result, stderr=log)
    transfer = {"exit_code": status.returncode, "archive": str(result_path),
                "archive_bytes": result_path.stat().st_size,
                "archive_sha256": sha256(result_path), "log": str(log_path),
                "log_sha256": sha256(log_path)}
    (output / "transfer.json").write_text(json.dumps(transfer, indent=2, sort_keys=True) + "\n")
    print(json.dumps(transfer, sort_keys=True), flush=True)
    raise SystemExit(status.returncode)


if __name__ == "__main__":
    main()
