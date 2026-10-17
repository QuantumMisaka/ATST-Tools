#!/usr/bin/env python3
"""Audit retained raw SAI acceptance archive and emit a compact case report."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import tarfile

REQUIRED_CASES = (
    "relax", "md_initial", "md_restart", "md_native", "neb_serial",
    "neb_parallel", "cp_baseline", "cp", "vibration", "lts_legacy",
    "sella", "ccqn", "irc",
)
REQUIRED_STAGES = ("venv", "install", "pip_check", "mpi_probe", "provenance") + REQUIRED_CASES


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def audit(directory: Path) -> dict:
    identity = json.loads((directory / "source-identity.json").read_text())
    transfer = json.loads((directory / "transfer.json").read_text())
    archive_path = Path(transfer["archive"])
    errors: list[str] = []
    if sha256(archive_path) != transfer["archive_sha256"]:
        errors.append("result archive SHA256 differs from transfer record")
    if sha256(Path(transfer["log"])) != transfer["log_sha256"]:
        errors.append("SAI log SHA256 differs from transfer record")
    if transfer["exit_code"] != 0:
        errors.append(f"streamed SAI process exit code {transfer['exit_code']}")
    cases = {}
    with tarfile.open(archive_path, "r:gz") as archive:
        members = {member.name: member for member in archive.getmembers()}

        def read(name: str) -> bytes:
            member = members.get(name)
            if member is None or not member.isfile():
                raise ValueError(f"missing archive member: {name}")
            with archive.extractfile(member) as stream:
                return stream.read()

        exit_code = read("exit-code.txt").decode().strip()
        if exit_code != "0":
            errors.append(f"node runtime exit code {exit_code}")
        embedded = json.loads(read("source-identity.json"))
        if embedded != identity:
            errors.append("embedded source identity differs from local input")
        stages = {}
        for line in read("stages.tsv").decode().splitlines():
            fields = line.split("\t")
            if len(fields) != 2 or fields[0] in stages:
                errors.append(f"malformed or duplicate stage: {line!r}")
                continue
            stages[fields[0]] = fields[1]
        for stage in REQUIRED_STAGES:
            if stages.get(stage) != "0":
                errors.append(f"stage {stage}: {stages.get(stage, 'missing')}")
        provenance = json.loads(read("harness/provenance.json"))
        if provenance["backend_source"] != "external":
            errors.append("installed backend is not external")
        if provenance["abacus_source_sha"] != identity["engine"]["source_sha"]:
            errors.append("engine source SHA differs from payload identity")
        if provenance["abacus_version_probe"]["returncode"] != 0:
            errors.append("engine version probe failed")
        for name, asset in identity["assets"].items():
            if provenance["assets"].get(name) != asset["sha256"]:
                errors.append(f"runtime asset differs: {name}")
        for key, identity_key in (("ATST_WHEEL", "atst"), ("PLUGIN_WHEEL", "abacuslite")):
            if provenance["wheels"][key]["sha256"] != identity[identity_key]["sha256"]:
                errors.append(f"installed {identity_key} wheel provenance differs")
            if provenance["wheels"][key]["matched_installed_python_files"] <= 0:
                errors.append(f"installed {identity_key} Python files were not byte-checked")
        for case in REQUIRED_CASES:
            root = "md_ase" if case in {"md_initial", "md_restart"} else case
            filename = ("evidence_initial.json" if case == "md_initial" else
                        "evidence_restart.json" if case == "md_restart" else "evidence.json")
            path = f"harness/cases/{root}/{filename}"
            if path not in members:
                errors.append(f"{case}: missing evidence")
                cases[case] = {"status": "missing"}
                continue
            evidence = json.loads(read(path))
            cases[case] = {"status": evidence["status"],
                           "api_status": evidence.get("checks", {}).get("api_status"),
                           "scientific_claim": evidence.get("scientific_claim"),
                           "error": evidence.get("error")}
            if evidence["status"] != "pass" or cases[case]["api_status"] != "complete":
                errors.append(f"{case}: {evidence['status']} / {cases[case]['api_status']}")
            if case in {"sella", "ccqn", "irc"} and evidence["status"] == "pass":
                stages_for_case = evidence["checks"]["workflow_manifest"]["stages"]
                cases[case]["optimizer_stages"] = stages_for_case
    return {"schema": "atst-runtime-acceptance-audit-v1", "archive": str(archive_path),
            "archive_sha256": transfer["archive_sha256"], "source_identity": identity,
            "stages": stages, "cases": cases, "errors": errors, "pass": not errors,
            "convergence_scope": "optimizer stage signals only; no scientific conclusion"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    result = audit(args.output_dir.resolve(strict=True))
    path = args.output_dir / "audit.json"
    path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"pass": result["pass"], "errors": result["errors"],
                      "audit": str(path)}, sort_keys=True))
    raise SystemExit(not result["pass"])


if __name__ == "__main__":
    main()
