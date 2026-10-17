#!/usr/bin/env python3
"""Real-ABACUS acceptance for installed external-only ATST and abacuslite.

This is an interface smoke, not a convergence or scientific-reference test.
Run from outside both source checkouts.  All generated files stay beside this
script.  The parallel case must be entered under an *outer* two-rank mpirun.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import hashlib
import importlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import traceback

HERE = Path(__file__).resolve().parent
CASES = HERE / "cases"
ASSETS = HERE / "assets"
FORBIDDEN = "atst_tools.external.ASE_interface.abacuslite"
TOP_LEVEL = "abacuslite"
EXPECTED_BACKEND = "external"
ATST_SHA = os.environ["ATST_SOURCE_SHA"]
PLUGIN_SHA = os.environ["PLUGIN_SOURCE_SHA"]


@contextmanager
def chdir(path):
    previous = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


import numpy as np  # noqa: E402
from ase.io import read, write  # noqa: E402
from atst_tools.api import RunOptions, run_workflow, validate_config  # noqa: E402
from atst_tools.calculators.abacuslite_backend import BACKEND_SOURCE, import_latestio, import_legacyio  # noqa: E402
from atst_tools.utils.io import read_structure  # noqa: E402
import atst_tools  # noqa: E402
import abacuslite  # noqa: E402


def require_backend():
    assert BACKEND_SOURCE == "external", BACKEND_SOURCE
    assert Path(atst_tools.__file__).resolve().is_relative_to(Path(sys.prefix).resolve())
    assert Path(abacuslite.__file__).resolve().is_relative_to(Path(sys.prefix).resolve())
    assert TOP_LEVEL in importlib.metadata.packages_distributions()
    legacy = importlib.import_module(FORBIDDEN + ".core")
    assert legacy.Abacus is abacuslite.Abacus, "historical import is not a thin forwarder"


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def dump(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name("." + path.name + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    tmp.replace(path)


def finite(value, label):
    arr = np.asarray(value, dtype=float)
    assert arr.size and np.isfinite(arr).all(), f"{label}: nonfinite/empty {arr.shape}"
    return arr


def atoms_facts(atoms, label, *, with_result=True):
    assert atoms is not None, f"{label}: missing atoms"
    xyz = finite(atoms.positions, label + " positions")
    assert xyz.shape == (len(atoms), 3)
    facts = {"natoms": len(atoms), "symbols": atoms.get_chemical_symbols(),
             "positions_angstrom": xyz.tolist(), "cell_angstrom": finite(atoms.cell, label + " cell").tolist()}
    if with_result:
        energy = float(atoms.get_potential_energy())
        forces = finite(atoms.get_forces(), label + " forces")
        assert math.isfinite(energy) and forces.shape == (len(atoms), 3)
        facts.update(energy_ev=energy, forces_ev_per_angstrom=forces.tolist(),
                     force_max_ev_per_angstrom=float(np.linalg.norm(forces, axis=1).max()))
    return facts


def newest_raw_log(workdir):
    logs = list(Path(workdir).rglob("running_scf.log")) + list(Path(workdir).rglob("running_relax.log"))
    assert logs, f"no ABACUS raw running log in {workdir}"
    return max(logs, key=lambda p: p.stat().st_mtime_ns)


def independent_raw_facts(log, natoms):
    """Parse literal ABACUS energy/force rows independently of abacuslite."""
    lines = Path(log).read_text(errors="replace").splitlines()
    number = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
    energy_rows = []
    fermi_rows = []
    force_tables = []
    for i, line in enumerate(lines):
        match = re.match(r"\s*E_KohnSham\s+(" + number + r")\s+(" + number + r")", line)
        if match:
            energy_rows.append(float(match.group(2).replace("D", "E")))
        match = re.match(r"\s*E_Fermi\s+(" + number + r")\s+(" + number + r")", line)
        if match:
            fermi_rows.append(float(match.group(2).replace("D", "E")))
        if re.search(r"TOTAL-FORCE\s*\(eV\s*/Angstrom\)", line, re.I):
            rows = []
            for next_line in lines[i + 1:]:
                parts = next_line.split()
                if len(parts) == 4 and re.fullmatch(r"[A-Z][a-z]?\d+", parts[0]):
                    try:
                        rows.append([float(v.replace("D", "E")) for v in parts[1:]])
                    except ValueError:
                        pass
                    if len(rows) == natoms:
                        force_tables.append(rows)
                        break
                elif rows:
                    break
    assert energy_rows, f"no literal E_KohnSham row in {log}"
    assert force_tables, f"no literal TOTAL-FORCE table in {log}"
    facts = {"energy_ev": energy_rows[-1], "forces_ev_per_angstrom": force_tables[-1],
             "energy_rows": len(energy_rows), "force_tables": len(force_tables)}
    if fermi_rows:
        facts["fermi_ev"] = fermi_rows[-1]
    return facts


def independent_md_dump_last(path, natoms):
    lines = Path(path).read_text(errors="replace").splitlines()
    frame_count = sum(line.strip() == "LATTICE_VECTORS" for line in lines)
    atom_rows = []
    pattern = re.compile(r"^\s*\d+\s+([A-Z][a-z]?)\s+([+-]?[\d.]+(?:[Ee][+-]?\d+)?)\s+([+-]?[\d.]+(?:[Ee][+-]?\d+)?)\s+([+-]?[\d.]+(?:[Ee][+-]?\d+)?)")
    for line in lines:
        match = pattern.match(line)
        if match:
            atom_rows.append((match.group(1), [float(match.group(i)) for i in (2, 3, 4)]))
    assert frame_count >= 2 and len(atom_rows) == frame_count * natoms
    final_rows = atom_rows[-natoms:]
    lattice_rows = [i for i, line in enumerate(lines) if line.strip() == "LATTICE_VECTORS"]
    last = lattice_rows[-1]
    alat_line = lines[last - 1].split()
    assert alat_line[0] == "LATTICE_CONSTANT:" and alat_line[-1] == "Angstrom"
    cell = np.asarray([lines[last + j].split() for j in (1, 2, 3)], float) * float(alat_line[1])
    return {"frame_count": frame_count, "cell_angstrom": cell.tolist(), "symbols": [r[0] for r in final_rows],
            "positions_angstrom": [r[1] for r in final_rows]}


def raw_compare(workdir, atoms, *, legacy=False, energy_key="energy_ev"):
    log = newest_raw_log(workdir)
    workflow = atoms_facts(atoms, "workflow")
    literal = independent_raw_facts(log, len(atoms))
    literal_energy_delta = abs(literal["energy_ev"] - workflow[energy_key])
    literal_force_delta = float(np.max(np.abs(np.asarray(literal["forces_ev_per_angstrom"]) -
                                             np.asarray(workflow["forces_ev_per_angstrom"]))))
    assert literal_energy_delta < 1e-5 and literal_force_delta < 1e-4
    eig_present = (log.parent / "eig_occ.txt").is_file()
    # The old vendored latestio unlinks eig_occ.txt on the workflow's first
    # read. Preserve this actual artifact difference; never reconstruct it.
    if not legacy and not eig_present:
        assert BACKEND_SOURCE == "vendored", "independent backend lost eig_occ.txt"
        parsed = None
        raw = {"energy_ev": literal["energy_ev"],
               "forces_ev_per_angstrom": literal["forces_ev_per_angstrom"]}
        reader_name = "literal_after_vendored_aux_removal"
    else:
        reader = import_legacyio() if legacy else import_latestio()
        parsed = reader.read_abacus_out(log, index=-1)
        if isinstance(parsed, list):
            parsed = parsed[-1]
        raw = atoms_facts(parsed, "raw ABACUS")
        assert raw["symbols"] == workflow["symbols"]
        reader_name = "legacyio" if legacy else "latestio"
    energy_diff = abs(raw["energy_ev"] - workflow[energy_key])
    force_diff = float(np.max(np.abs(np.asarray(raw["forces_ev_per_angstrom"]) -
                                     np.asarray(workflow["forces_ev_per_angstrom"]))))
    assert energy_diff < 1e-5, f"workflow/raw energy drift {energy_diff} eV"
    assert force_diff < 1e-4, f"workflow/raw force drift {force_diff} eV/Å"
    fermi = parsed.calc.get_fermi_level() if parsed is not None else None
    if fermi is None:
        fermi = literal.get("fermi_ev")
    return {"reader": reader_name, "eig_occ_present_after_workflow": eig_present,
            "log": str(log.relative_to(HERE)),
            "log_sha256": sha256(log), "energy_abs_delta_ev": energy_diff,
            "force_max_abs_delta_ev_per_angstrom": force_diff,
            "literal_raw": literal, "literal_energy_abs_delta_ev": literal_energy_delta,
            "literal_force_max_abs_delta_ev_per_angstrom": literal_force_delta,
            "raw_energy_ev": raw["energy_ev"],
            "raw_fermi_ev": None if fermi is None else float(fermi)}


def provenance():
    require_backend()
    wheels = {}
    keys = ("ATST_WHEEL", "PLUGIN_WHEEL")
    for key in keys:
        value = os.environ.get(key)
        if not value:
            raise RuntimeError(f"{key} is required for wheel provenance")
        path = Path(value).resolve(strict=True)
        wheels[key] = {"path": str(path), "sha256": sha256(path)}
    assert wheels["ATST_WHEEL"]["sha256"] == os.environ["EXPECTED_ATST_SHA256"]
    assert wheels["PLUGIN_WHEEL"]["sha256"] == os.environ["EXPECTED_PLUGIN_SHA256"]
    import zipfile
    packages = (("ATST_WHEEL", "atst_tools"), ("PLUGIN_WHEEL", "abacuslite"))
    for key, package in packages:
        installed = Path(importlib.import_module(package).__file__).resolve().parent
        assert installed.is_relative_to(Path(sys.prefix).resolve()), installed
        with zipfile.ZipFile(wheels[key]["path"]) as archive:
            names = [n for n in archive.namelist() if n.startswith(package + "/") and n.endswith(".py")]
            for name in names:
                assert (installed.parent / name).read_bytes() == archive.read(name), name
        wheels[key]["matched_installed_python_files"] = len(names)
    command = os.environ.get("ABACUS_COMMAND")
    if not command:
        raise RuntimeError("ABACUS_COMMAND is required")
    version_cmd = os.environ.get("ABACUS_VERSION_COMMAND")
    version = None
    if version_cmd:
        output = subprocess.run(shlex.split(version_cmd), capture_output=True, text=True, timeout=30)
        version = {"command": version_cmd, "returncode": output.returncode,
                   "stdout": output.stdout, "stderr": output.stderr}
        assert output.returncode == 0, "ABACUS version probe failed"
    return {"utc": datetime.now(timezone.utc).isoformat(), "python": sys.version,
            "atst": {"version": importlib.metadata.version("atst-tools"), "path": atst_tools.__file__, "source_sha": ATST_SHA},
            "abacuslite": {"version": importlib.metadata.version("abacuslite"), "path": abacuslite.__file__, "source_sha": PLUGIN_SHA},
            "backend_source": BACKEND_SOURCE, "wheels": wheels,
            "abacus_command": command, "abacus_source_sha": os.environ.get("ABACUS_SOURCE_SHA"),
            "abacus_version_probe": version, "assets": {p.name: sha256(p) for p in ASSETS.iterdir() if p.is_file()}}


def abacus_config(case, *, species="H", command=None, cp=False):
    command = command or os.environ.get("ABACUS_COMMAND", "abacus")
    prefix = species.lower()
    params = {"calculation": "scf", "basis_type": "lcao", "ks_solver": "cusolver",
              "dft_functional": "pbe", "ecutwfc": 100, "scf_thr": 1e-7,
              "scf_nmax": 100, "cal_force": 1, "cal_stress": 0,
              "pseudo_dir": str(ASSETS), "orbital_dir": str(ASSETS),
              "pseudopotentials": {species: f"{prefix}.upf"},
              "basissets": {species: f"{prefix}.orb"}}
    if cp:
        params.update({"suffix": "CP", "symmetry": 0, "nspin": 1, "nbands": 130,
                       "scf_thr": 1e-9, "scf_nmax": 200, "smearing_method": "gauss",
                       "smearing_sigma": 0.015, "mixing_type": "broyden", "mixing_beta": 0.1,
                       "out_pot": 2, "out_chg": "1 12", "gate_flag": 1, "efield_dir": 0,
                       "zgate": 0.7, "block_down": 0.65, "block_up": 0.75,
                       "block_height": 0.01, "block": 1, "nelec": 217,
                       "efield_flag": 1, "dip_cor_flag": 1, "efield_pos_max": 0.8,
                       "efield_pos_dec": 0.1, "efield_amp": 0})
    return {"name": "abacus", "abacus": {"command": command, "mpi": 1, "omp": 1,
            "directory": "abacus_run", "kpts": [1, 1, 1], "parameters": params}}


def config_for(case, target_mu=None):
    if case in ("relax", "lts_legacy"):
        calc = {"type": "relax", "init_structure": str(ASSETS / "h2.stru"),
                "fmax": 0.02, "max_steps": 1, "optimizer": "FIRE",
                "trajectory": "relax.traj", "logfile": "relax.log"}
        command = os.environ.get("ABACUS_LTS_COMMAND") if case == "lts_legacy" else None
        if case == "lts_legacy" and not command:
            raise RuntimeError("ABACUS_LTS_COMMAND is required for LTS 3.10.1 reader check")
        return {"calculation": calc, "calculator": abacus_config(case, command=command)}
    if case in ("md_initial", "md_restart"):
        return {"calculation": {"type": "md", "driver": "ase", "init_structure": str(ASSETS / "h2.stru"),
                "ensemble": "nvt", "algorithm": "langevin", "steps": 1,
                "timestep_fs": 0.5, "temperature_K": 300.0, "friction_fs_inv": 0.01,
                "seed": 1729, "trajectory": "md.traj", "logfile": "md.log", "summary_file": "md_summary.json",
                "final_structure": "md_final.traj"}, "calculator": abacus_config(case)}
    if case == "md_native":
        config = abacus_config(case)
        config["abacus"]["parameters"].update({"calculation": "md", "symmetry": 0,
                                                   "md_type": "nvt", "md_nstep": 2,
                                                   "md_dt": 0.5, "md_tfirst": 300,
                                                   "md_seed": 1729})
        return {"calculation": {"type": "md", "driver": "abacus_native",
                "init_structure": str(ASSETS / "h2.stru"), "steps": 2,
                "directory": "native_run", "trajectory": "md_native.traj",
                "summary_file": "md_native_summary.json", "final_structure": "md_native_final.traj",
                "poll_interval_seconds": 1.0, "timeout_seconds": 900},
                "calculator": config}
    if case in ("neb_serial", "neb_parallel"):
        return {"calculation": {"type": "neb", "init_chain": str(ASSETS / "h2_chain.traj"),
                "climb": False, "two_stage": False, "fmax": 0.05, "k": 0.1,
                "algorism": "improvedtangent", "parallel": case == "neb_parallel",
                "endpoint_singlepoint": "always", "optimizer": "FIRE", "max_steps": 1},
                "calculator": abacus_config(case)}
    if case == "cp_baseline":
        return {"calculation": {"type": "relax", "init_structure": str(ASSETS / "pt.stru"),
                "fmax": 1e6, "max_steps": 1, "optimizer": "FIRE",
                "trajectory": "relax.traj", "logfile": "relax.log"},
                "calculator": abacus_config(case, species="Pt", cp=True)}
    if case == "cp":
        assert target_mu is not None
        return {"calculation": {"type": "constant_potential", "init_structure": str(ASSETS / "pt.stru"),
                "directory": "cp_run", "results_file": "cp_results.json", "log_file": "cp.log",
                "checkpoint_file": "cp_checkpoint.json"},
                "calculator": {**abacus_config(case, species="Pt", cp=True),
                    "constant_potential": {"energy_boundary": "compensated_gate",
                    "target_mu_ev": float(target_mu), "reference_electrode": "custom",
                    "reference_electrons": 216, "potential_tolerance_v": 0.02,
                    "max_iterations": 2, "nelec_min": 216.5, "nelec_max": 217.5,
                    "nelec_step_max": 0.2, "vacuum_axis": 0}}}
    if case == "vibration":
        return {"calculation": {"type": "vibration", "init_structure": str(ASSETS / "h2.stru"),
                "indices": [0], "delta": 0.01, "nfree": 2, "name": "vib",
                "results_file": "vibration_results.json", "validation_file": "ts_validation.json",
                "thermochemistry": {"model": "harmonic", "temperature": 300.0}},
                "calculator": abacus_config(case)}
    if case == "sella":
        return {"calculation": {"type": "sella", "init_structure": str(ASSETS / "h2.stru"),
                "trajectory": "sella.traj", "fmax": 0.02, "max_steps": 1,
                "eta": 0.002, "order": 1, "record_events": True},
                "calculator": abacus_config(case)}
    if case == "ccqn":
        return {"calculation": {"type": "ccqn", "init_structure": str(ASSETS / "h2.stru"),
                "trajectory": "ccqn.traj", "logfile": "ccqn.log",
                "final_structure": "ccqn_final.extxyz", "artifact_manifest": "atst_artifacts.json",
                "fmax": 0.02, "max_steps": 1, "e_vector_method": "ic",
                "reactive_bonds": "1-2", "accept_initial_converged": False},
                "calculator": abacus_config(case)}
    if case == "irc":
        return {"calculation": {"type": "irc", "backend": "descent",
                "init_structure": str(ASSETS / "h2.stru"), "mode_vector": str(ASSETS / "h2_mode.npy"),
                "trajectory": "irc.traj", "normalized_trajectory": "irc_normalized.traj",
                "direction": "both", "descent_delta": 0.05, "fmax": 0.02,
                "max_steps": 1, "artifact_manifest": "atst_artifacts.json"},
                "calculator": abacus_config(case)}
    raise ValueError(case)


def prepare(source):
    require_backend()
    source = Path(source).resolve(strict=True)
    if not (source / "data").is_dir():
        raise RuntimeError("--source-examples must be ATST examples directory")
    ASSETS.mkdir(exist_ok=True)
    mappings = {
        "h2.stru": source / "11_vibration_ideal_gas_H2/inputs/h2.stru",
        "pt.stru": source / "19_constant_potential_Pt/inputs/STRU",
        "h.upf": source / "data/H_ONCV_PBE-1.0.upf",
        "h.orb": source / "data/H_gga_6au_100Ry_2s1p.orb",
        "pt.upf": source / "data/Pt_ONCV_PBE-1.0.upf",
        "pt.orb": source / "data/Pt_gga_7au_100Ry_4s2p2d1f.orb",
    }
    for name, src in mappings.items():
        dst = ASSETS / name
        if dst.exists() and sha256(dst) != sha256(src):
            raise RuntimeError(f"existing asset differs: {dst}")
        shutil.copy2(src, dst)
    h2 = read_structure(str(ASSETS / "h2.stru"))
    assert len(h2) == 2 and h2.get_chemical_symbols() == ["H", "H"]
    chain = []
    for delta in (0.0, 0.02, 0.04, 0.06):
        image = h2.copy()
        image.positions[1, 0] += delta
        chain.append(image)
    write(str(ASSETS / "h2_chain.traj"), chain)
    mode = np.zeros((2, 3), dtype=float)
    mode[0, 0], mode[1, 0] = -1.0, 1.0
    np.save(ASSETS / "h2_mode.npy", mode)
    pt = read_structure(str(ASSETS / "pt.stru"))
    assert len(pt) == 12 and set(pt.get_chemical_symbols()) == {"Pt"}
    # Check the real public schema before cluster execution; CP target is
    # replaced with a same-engine value after cp_baseline completes.
    checked = {}
    for case in ("relax", "md_initial", "md_restart", "md_native", "neb_serial", "neb_parallel", "cp_baseline", "cp", "vibration", "sella", "ccqn", "irc"):
        config = config_for(case, target_mu=0.0)
        validate_config(config)
        checked[case] = config["calculation"]["type"]
    dump(HERE / "prepared.json", {"schema": "atst-real-smoke-prepared-v1", "source_examples": str(source),
                                  "source_sha": ATST_SHA, "assets": {p.name: sha256(p) for p in ASSETS.iterdir()},
                                  "validated_config_types": checked,
                                  "note": "H2 bond interpolation is an interface smoke, not a reaction path."})
    print(json.dumps(checked, sort_keys=True))


def case_root(case):
    return CASES / ("md_ase" if case in ("md_initial", "md_restart") else case)


def snapshot_md_initial(workdir, raw_log):
    """Freeze initial MD output before the restart reuses its work directory."""
    snapshot = workdir / "md_initial_snapshot"
    snapshot.mkdir(exist_ok=False)
    source_paths = {
        "raw_log": raw_log,
        "trajectory": workdir / "md.traj",
        "summary": workdir / "md_summary.json",
        "final_structure": workdir / "md_final.traj",
        "input": raw_log.parent.parent / "INPUT",
        "structure": raw_log.parent.parent / "STRU",
    }
    recorded = {}
    for role, source in source_paths.items():
        source = Path(source).resolve(strict=True)
        target = snapshot / role / source.name
        target.parent.mkdir(parents=True, exist_ok=False)
        original_sha = sha256(source)
        shutil.copy2(source, target)
        assert sha256(target) == original_sha, f"MD initial snapshot changed: {role}"
        recorded[role] = {"path": str(target.relative_to(HERE)), "sha256": original_sha}
    return recorded


def verify_md_initial_snapshot(prior):
    recorded = prior["checks"]["initial_snapshot"]
    for role in ("raw_log", "trajectory", "summary", "final_structure", "input", "structure"):
        entry = recorded[role]
        path = HERE / entry["path"]
        assert path.is_file() and sha256(path) == entry["sha256"], f"MD initial snapshot lost: {role}"
    assert recorded["raw_log"]["path"] == prior["checks"]["raw_log"]
    return recorded


def mpi_world(case):
    if case != "neb_parallel":
        return None, True
    from mpi4py import MPI
    from ase.parallel import world
    assert MPI.COMM_WORLD.size == world.size, "ASE and site MPI disagree"
    assert world.size == 2, f"image-parallel NEB needs exactly 2 outer ranks, got {world.size}"
    return world, world.rank == 0


def run_case(case):
    require_backend()
    world, root = mpi_world(case)
    workdir = case_root(case)
    workdir.mkdir(parents=True, exist_ok=True)
    if case == "md_restart" and not (workdir / "evidence_initial.json").is_file():
        raise RuntimeError("run md_initial before md_restart")
    if case == "md_restart":
        verify_md_initial_snapshot(json.loads((workdir / "evidence_initial.json").read_text()))
    target = None
    if case == "cp":
        baseline = json.loads((case_root("cp_baseline") / "evidence.json").read_text())
        target = baseline["checks"]["derived_mu_ev"]
    config = config_for(case, target)
    validate_config(config)
    if root:
        dump(workdir / ("config_" + case + ".json"), config)
    if world:
        world.barrier()
    evidence = {"schema": "atst-real-smoke-case-v1", "case": case,
                "started_utc": datetime.now(timezone.utc).isoformat(),
                "scientific_claim": "bounded interface smoke only", "status": "running"}
    try:
        if case in ("md_initial", "md_restart"):
            # ATST seeds initial velocities; ASE Langevin uses global NumPy RNG.
            np.random.seed(1729)
        with chdir(workdir):
            result = run_workflow(config, RunOptions(restart=case == "md_restart", world=world, profiles=case.startswith("neb")))
            checks = {"workflow": result.workflow, "api_status": result.status,
                      "is_root": result.is_root, "manifest": result.artifact_manifest,
                      "backend_source": BACKEND_SOURCE}
            if root:
                assert result.status == "complete", result.status
                document = result.to_document(workdir)
                dump(workdir / ("api_" + case + ".json"), document)
                if case in ("relax", "lts_legacy", "cp_baseline"):
                    final = read("final_relaxed.traj")
                    checks["final"] = atoms_facts(final, "final relax")
                    checks["raw_compare"] = raw_compare(workdir, final, legacy=case == "lts_legacy")
                    if case == "cp_baseline":
                        from atst_tools.utils.gate_compensation import read_gate_compensation
                        raw_log = newest_raw_log(workdir)
                        evaldir = raw_log.parent.parent
                        compensation = read_gate_compensation(evaldir, 217.0, 216.0)
                        checks["compensation"] = compensation
                        assert checks["raw_compare"]["raw_fermi_ev"] is not None
                        checks["derived_mu_ev"] = float(checks["raw_compare"]["raw_fermi_ev"] +
                                                        compensation["compensation_derivative_ev"])
                        finite(checks["derived_mu_ev"], "derived CP mu")
                elif case in ("md_initial", "md_restart"):
                    frames = read("md.traj", index=":")
                    assert len(frames) >= 2
                    checks["frame_count"] = len(frames)
                    checks["last_frame"] = atoms_facts(frames[-1], "MD last", with_result=False)
                    checks["summary"] = json.loads(Path("md_summary.json").read_text())
                    assert checks["summary"]["steps"] == config["calculation"]["steps"]
                    assert len(frames) == checks["summary"]["steps"] + 1
                    finite(checks["summary"]["final_energy_eV"], "MD summary energy")
                    finite(checks["summary"].get("final_temperature_K", 0.0), "MD summary temperature")
                    if case == "md_restart":
                        prior = json.loads((workdir / "evidence_initial.json").read_text())
                        previous_last = np.asarray(prior["checks"]["last_frame"]["positions_angstrom"])
                        # ASE's trajectory writer may replace the file on a
                        # new run; continuity is the first restart frame.
                        assert np.allclose(frames[0].positions, previous_last, atol=1e-12)
                        previous_energy = float(prior["checks"]["summary"]["final_energy_eV"])
                        boundary_energy = float(frames[0].get_potential_energy())
                        assert abs(boundary_energy - previous_energy) < 1e-5
                        checks["restart_boundary_energy_abs_delta_ev"] = abs(boundary_energy - previous_energy)
                        checks["restart_boundary_max_delta_angstrom"] = float(
                            np.max(np.abs(frames[0].positions - previous_last)))
                        assert np.max(np.abs(frames[-1].positions - previous_last)) > 0, "MD restart made no geometrical progress"
                    raw_log = newest_raw_log(workdir)
                    literal = independent_raw_facts(raw_log, 2)
                    assert abs(literal["energy_ev"] - float(checks["summary"]["final_energy_eV"])) < 1e-5
                    checks["literal_raw_last"] = literal
                    if case == "md_initial":
                        checks["initial_snapshot"] = snapshot_md_initial(workdir, raw_log)
                        checks["raw_log"] = checks["initial_snapshot"]["raw_log"]["path"]
                    else:
                        checks["raw_log"] = str(raw_log.relative_to(HERE))
                    checks["raw_log_sha256"] = sha256(HERE / checks["raw_log"])
                elif case == "md_native":
                    summary = json.loads(Path("md_native_summary.json").read_text())
                    assert summary["driver"] == "abacus_native" and summary["returncode"] == 0
                    frames = read("md_native.traj", index=":")
                    assert len(frames) >= 2 and summary["frames"] == len(frames)
                    checks["summary"] = summary
                    checks["frame_count"] = len(frames)
                    checks["last_frame"] = atoms_facts(frames[-1], "native MD last", with_result=False)
                    dump_path = Path("native_run/OUT.ABACUS/MD_dump")
                    assert dump_path.is_file() and dump_path.stat().st_size > 0
                    from atst_tools.calculators.abacuslite_backend import import_legacyio
                    raw_frames = import_legacyio().read_traj_from_md_dump(dump_path)
                    assert len(raw_frames) == len(frames)
                    raw_last = finite(raw_frames[-1]["coords"], "MD_dump last positions")
                    delta = float(np.max(np.abs(raw_last - frames[-1].positions)))
                    assert delta < 1e-5
                    literal_dump = independent_md_dump_last(dump_path, 2)
                    assert literal_dump["frame_count"] == len(frames)
                    assert literal_dump["symbols"] == ["H", "H"]
                    literal_delta = float(np.max(np.abs(np.asarray(literal_dump["positions_angstrom"]) - frames[-1].positions)))
                    assert literal_delta < 1e-5
                    cell_delta = float(np.max(np.abs(np.asarray(literal_dump["cell_angstrom"]) - frames[-1].cell)))
                    assert cell_delta < 1e-5
                    assert all(frame.pbc.all() for frame in frames)
                    checks["literal_cell_max_delta_angstrom"] = cell_delta
                    checks["periodic_frames"] = True
                    checks["md_dump_sha256"] = sha256(dump_path)
                    checks["dump_last_position_max_delta_angstrom"] = delta
                    checks["literal_md_dump_last"] = literal_dump
                    checks["literal_dump_last_position_max_delta_angstrom"] = literal_delta
                elif case.startswith("neb_"):
                    assert result.final_images is not None and len(result.final_images) == 4
                    checks["images"] = [atoms_facts(img, f"NEB image {i}", with_result=False)
                                        for i, img in enumerate(result.final_images)]
                    assert len(result.profiles) == 4
                    checks["profiles"] = list(result.profiles)
                    # The optional API profiles may be null for an image owned
                    # by another MPI rank. The persisted trajectory freezes the
                    # globally gathered band energy/forces and is authoritative.
                    stored = read("neb.traj", index=":")
                    assert len(stored) >= 4 and len(stored) % 4 == 0
                    final_band = stored[-4:]
                    checks["persisted_band"] = [atoms_facts(frame, f"persisted NEB {i}")
                                                for i, frame in enumerate(final_band)]
                    for i, frame in enumerate(final_band):
                        assert np.allclose(frame.positions, result.final_images[i].positions, atol=1e-10)
                        profile = checks["profiles"][i]
                        if profile["energy_eV"] is not None:
                            assert abs(profile["energy_eV"] - frame.get_potential_energy()) < 1e-5
                            assert np.allclose(profile["forces"], frame.get_forces(), atol=1e-5)
                        else:
                            assert case == "neb_parallel"
                    raw_outputs = []
                    for log in workdir.rglob("running_scf.log"):
                        source_stru = log.parent.parent / "STRU"
                        assert source_stru.is_file(), source_stru
                        input_atoms = read_structure(str(source_stru))
                        raw_facts = independent_raw_facts(log, 2)
                        if BACKEND_SOURCE == "external":
                            assert (log.parent / "eig_occ.txt").is_file(), log
                            parsed = import_latestio().read_abacus_out(log, index=-1)
                            assert np.allclose(parsed.positions, input_atoms.positions, atol=1e-5)
                            assert abs(parsed.get_potential_energy() - raw_facts["energy_ev"]) < 1e-5
                        raw_outputs.append((log, input_atoms, raw_facts))
                    checked_logs = []
                    for frame in final_band:
                        matches = [(log, facts) for log, parsed, facts in raw_outputs
                                   if np.allclose(parsed.positions, frame.positions, atol=1e-5)
                                   and np.allclose(parsed.cell, frame.cell, atol=1e-5)]
                        assert matches, "no raw ABACUS output matches a final band image"
                        assert any(abs(facts["energy_ev"] - frame.get_potential_energy()) < 1e-5
                                   and np.allclose(facts["forces_ev_per_angstrom"], frame.get_forces(), atol=1e-4)
                                   for _, facts in matches)
                        checked_logs.append([str(log.relative_to(HERE)) for log, _ in matches])
                    checks["persisted_band_raw_matches"] = checked_logs
                    checks["raw_logs"] = [str(p.relative_to(HERE)) for p in workdir.rglob("running_scf.log")]
                    assert len(checks["raw_logs"]) >= 2, "missing image ABACUS logs"
                    checks["image_count"] = 4
                    checks["mpi_world_size"] = world.size if world else 1
                elif case == "cp":
                    payload = json.loads(Path("cp_results.json").read_text())
                    assert payload["energy_boundary"] == "compensated_gate"
                    assert len(payload["points"]) == 1
                    point = payload["points"][0]
                    assert point["cp_converged"] is True
                    for key in ("energy", "raw_energy", "mu_calc", "residual_mu_eV", "nelec"):
                        finite(point[key], "CP " + key)
                    assert abs(point["residual_mu_eV"]) <= 0.02
                    forces = finite(point["forces"], "CP forces")
                    assert forces.shape == (12, 3)
                    target_mu = float(config["calculator"]["constant_potential"]["target_mu_ev"])
                    expected_omega = float(point["raw_energy"]) - target_mu * (float(point["nelec"]) - 216.0)
                    assert abs(float(point["energy"]) - expected_omega) < 1e-6
                    checks["cp_point"] = point
                    checks["omega_bookkeeping_abs_delta_ev"] = abs(float(point["energy"]) - expected_omega)
                    checks["target_from_baseline_ev"] = target_mu
                    raw_log = newest_raw_log(workdir)
                    literal = independent_raw_facts(raw_log, 12)
                    eig_present = (raw_log.parent / "eig_occ.txt").is_file()
                    if eig_present:
                        raw = import_latestio().read_abacus_out(raw_log, index=-1)
                        if isinstance(raw, list):
                            raw = raw[-1]
                        raw_facts = atoms_facts(raw, "CP raw ABACUS")
                    else:
                        assert BACKEND_SOURCE == "vendored", "independent CP lost eig_occ.txt"
                        raw_facts = {"energy_ev": literal["energy_ev"],
                                     "forces_ev_per_angstrom": literal["forces_ev_per_angstrom"]}
                    assert abs(raw_facts["energy_ev"] - float(point["raw_energy"])) < 1e-5
                    assert abs(literal["energy_ev"] - float(point["raw_energy"])) < 1e-5
                    force_delta = float(np.max(np.abs(np.asarray(raw_facts["forces_ev_per_angstrom"]) - forces)))
                    literal_force_delta = float(np.max(np.abs(np.asarray(literal["forces_ev_per_angstrom"]) - forces)))
                    assert force_delta < 1e-4
                    assert literal_force_delta < 1e-4
                    checks["raw_compare"] = {"log": str(raw_log.relative_to(HERE)), "log_sha256": sha256(raw_log),
                                             "eig_occ_present_after_workflow": eig_present,
                                             "raw_energy_abs_delta_ev": abs(raw_facts["energy_ev"] - float(point["raw_energy"])),
                                             "raw_force_max_abs_delta_ev_per_angstrom": force_delta,
                                             "literal_raw": literal,
                                             "literal_energy_abs_delta_ev": abs(literal["energy_ev"] - float(point["raw_energy"])),
                                             "literal_force_max_abs_delta_ev_per_angstrom": literal_force_delta}
                elif case == "vibration":
                    vib = json.loads(Path("vibration_results.json").read_text())
                    frequencies = finite(vib["frequencies"], "vibration frequencies")
                    imaginary = finite(vib["imaginary_frequencies"], "vibration imaginary frequencies")
                    assert frequencies.shape == (3,) and imaginary.shape == (3,)
                    assert vib["indices"] == [0] and vib["n_atoms"] == 2
                    finite(vib["zpe"], "vibration zpe")
                    finite(vib["thermo"]["free_energy"], "vibration free energy")
                    cache = sorted(Path("vib").glob("cache.*.json"))
                    assert len(cache) >= 6, f"missing displacement cache: {len(cache)}"
                    cache_forces = []
                    from ase.io.jsonio import read_json
                    for entry in cache:
                        payload = read_json(entry, always_array=False)
                        rows = finite(payload["forces"], entry.name + " forces")
                        assert rows.shape == (2, 3)
                        cache_forces.append(rows.tolist())
                    raw_log = newest_raw_log(workdir)
                    literal = independent_raw_facts(raw_log, 2)
                    assert any(np.allclose(rows, literal["forces_ev_per_angstrom"], atol=1e-4)
                               for rows in cache_forces), "final raw force absent from vibration cache"
                    checks["vibration"] = vib
                    checks["displacement_cache_count"] = len(cache)
                    checks["raw_last"] = {"log": str(raw_log.relative_to(HERE)),
                                          "log_sha256": sha256(raw_log), "literal_raw": literal}
                elif case in ("sella", "ccqn", "irc"):
                    manifest_path = Path("atst_artifacts.json")
                    assert manifest_path.is_file(), f"{case}: missing artifact manifest"
                    manifest = json.loads(manifest_path.read_text())
                    assert manifest["schema_version"] == "atst-artifacts-v1"
                    assert manifest["workflow"] == case
                    assert manifest["stages"] and all(
                        stage["status"] == "complete" for stage in manifest["stages"]
                    ), manifest["stages"]
                    for artifact in manifest["artifacts"]:
                        assert Path(artifact["path"]).is_file(), artifact
                    checks["workflow_manifest"] = manifest
                    if case == "irc":
                        assert len(manifest["stages"]) == 2
                        assert {s["direction"] for s in manifest["stages"]} == {"forward", "reverse"}
                        frames = read("irc.traj", index=":")
                        assert len(frames) == 2
                        assert len(read("irc_normalized.traj", index=":")) == 2
                        assert not np.allclose(frames[0].positions, frames[1].positions, atol=1e-8)
                        checks["endpoint_frames"] = [atoms_facts(frame, f"IRC endpoint {i}")
                                                      for i, frame in enumerate(frames)]
                    else:
                        trajectory = "sella.traj" if case == "sella" else "ccqn.traj"
                        frames = read(trajectory, index=":")
                        assert len(frames) >= 1
                        checks["candidate_frame"] = atoms_facts(frames[-1], case + " candidate")
                        checks["frame_count"] = len(frames)
                        if case == "ccqn":
                            assert Path("ccqn_final.extxyz").is_file()
                            assert Path("ccqn_mode_manifest.json").is_file()
                            checks["mode_manifest"] = json.loads(Path("ccqn_mode_manifest.json").read_text())
                    raw_log = newest_raw_log(workdir)
                    checks["literal_raw_last"] = independent_raw_facts(raw_log, 2)
                    checks["literal_raw_log"] = str(raw_log.relative_to(HERE))
                    checks["literal_raw_log_sha256"] = sha256(raw_log)
                evidence["checks"] = checks
                evidence["status"] = "pass"
    except BaseException as exc:
        if root:
            evidence["status"] = "fail"
            evidence["error"] = {"type": type(exc).__name__, "message": str(exc),
                                 "traceback": traceback.format_exc()}
        raise
    finally:
        if root:
            evidence["finished_utc"] = datetime.now(timezone.utc).isoformat()
            filename = "evidence_initial.json" if case == "md_initial" else "evidence_restart.json" if case == "md_restart" else "evidence.json"
            dump(workdir / filename, evidence)
            print(json.dumps({"case": case, "status": evidence["status"], "evidence": str(workdir / filename)}), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subs = parser.add_subparsers(dest="action", required=True)
    prep = subs.add_parser("prepare")
    prep.add_argument("--source-examples", required=True)
    subs.add_parser("provenance")
    run = subs.add_parser("run-case")
    run.add_argument("case", choices=("relax", "md_initial", "md_restart", "md_native", "neb_serial", "neb_parallel", "cp_baseline", "cp", "lts_legacy", "vibration", "sella", "ccqn", "irc"))
    args = parser.parse_args()
    if args.action == "prepare":
        prepare(args.source_examples)
    elif args.action == "provenance":
        dump(HERE / "provenance.json", provenance())
        print(HERE / "provenance.json")
    else:
        if not (HERE / "provenance.json").is_file():
            raise RuntimeError("run provenance before any case")
        recorded = json.loads((HERE / "provenance.json").read_text())
        current = provenance()
        assert recorded["wheels"] == current["wheels"], "wheel identity changed since provenance"
        assert recorded["abacus_command"] == current["abacus_command"], "ABACUS_COMMAND changed since provenance"
        run_case(args.case)


if __name__ == "__main__":
    main()
