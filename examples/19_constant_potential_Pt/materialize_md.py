"""Materialize a fixed-cell Pt constant-potential MD input without running ABACUS.

The target chemical potential must be measured or chosen for this exact LCAO,
gate, dipole, and blocking-potential Hamiltonian. It is deliberately required
on the command line rather than copied from another fixture.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

from ase.io import write
from ruamel.yaml import YAML

from materialize import materialize
from atst_tools.utils.io import read_abacus_stru


def materialize_md(destination: Path, target_mu_ev: float, *, steps: int = 5) -> Path:
    """Create a constrained Pt slab, pinned assets, and a CP MD YAML config.

    Args:
        destination: New self-contained example directory.
        target_mu_ev: Explicit chemical-potential target for this Hamiltonian.
        steps: Number of MD integration steps.

    Returns:
        Path to the generated configuration.

    Raises:
        ValueError: If the target or step count is invalid.
        FileExistsError: If the destination already exists.
    """
    if not math.isfinite(target_mu_ev):
        raise ValueError("target_mu_ev must be finite")
    if steps < 1:
        raise ValueError("steps must be positive")
    destination = destination.resolve()
    destination.mkdir(parents=True, exist_ok=False)
    fixture = materialize(
        destination / "fixture", "lcao", 217.0, dipole=True, block_height=2.0
    )
    stru_path = fixture / "STRU"
    lines = stru_path.read_text(encoding="utf-8").splitlines()
    atom_lines = [index for index, line in enumerate(lines) if line.strip().endswith("1 1 1")]
    if len(atom_lines) != 12:
        raise ValueError("The pinned Pt fixture must contain 12 initially mobile atoms")
    for index in atom_lines[:8]:
        lines[index] = lines[index].rsplit("1 1 1", 1)[0] + "0 0 0"
    stru_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    atoms = read_abacus_stru(stru_path)
    write(fixture / "md_init.traj", atoms)

    # These are the settings of materialize.py's pinned INPUT. The CP factory
    # supplies out_chg='1 12' and a fresh nelec at each electronic evaluation.
    parameters = {
        "calculation": "scf",
        "suffix": "CP",
        "dft_functional": "pbe",
        "symmetry": 0,
        "cal_force": 1,
        "cal_stress": 0,
        "ecutwfc": 100,
        "scf_nmax": 200,
        "scf_thr": 1e-9,
        "basis_type": "lcao",
        "ks_solver": "cusolver",
        "nspin": 1,
        "nbands": 130,
        "smearing_method": "gauss",
        "smearing_sigma": 0.015,
        "mixing_type": "broyden",
        "mixing_beta": 0.1,
        "out_pot": 2,
        "out_chg": "1 12",
        "gate_flag": 1,
        "efield_dir": 0,
        "zgate": 0.7,
        "block_down": 0.65,
        "block_up": 0.75,
        "block_height": 2.0,
        "block": 1,
        "efield_flag": 1,
        "dip_cor_flag": 1,
        "efield_pos_max": 0.8,
        "efield_pos_dec": 0.1,
        "efield_amp": 0,
        "nelec": 217.0,
    }
    config = {
        "calculation": {
            "type": "md",
            "driver": "ase",
            "init_structure": "fixture/md_init.traj",
            "ensemble": "nvt",
            "algorithm": "langevin",
            "steps": steps,
            "timestep_fs": 0.25,
            "temperature_K": 100.0,
            "friction_fs_inv": 0.01,
            "seed": 19,
            "directory": "run_md",
            "trajectory": "run_md/md.traj",
            "logfile": "run_md/md.log",
            "summary_file": "run_md/md_summary.json",
            "final_structure": "run_md/md_final.traj",
            "artifact_manifest": "run_md/atst_artifacts.json",
        },
        "calculator": {
            "name": "abacus",
            "abacus": {
                "command": "abacus",
                "mpi": 1,
                "omp": 8,
                "directory": "run_md",
                "kpts": [1, 1, 1],
                "pseudo_dir": "fixture",
                "orbital_dir": "fixture",
                "pseudopotentials": {"Pt": "Pt_ONCV_PBE-1.0.upf"},
                "basissets": {"Pt": "Pt_gga_7au_100Ry_4s2p2d1f.orb"},
                "parameters": parameters,
            },
            "constant_potential": {
                "energy_boundary": "compensated_gate",
                "reference_electrode": "custom",
                "target_mu_ev": float(target_mu_ev),
                "reference_electrons": 216.0,
                "potential_tolerance_v": 1e-4,
                "max_iterations": 12,
                "capacitance_initial": 0.05,
                "capacitance_unit": "e/V",
                "nelec_min": 216.5,
                "nelec_max": 217.5,
                "nelec_step_max": 0.1,
                "vacuum_axis": 0,
            },
        },
    }
    path = destination / "config_md.yaml"
    with path.open("w", encoding="utf-8") as handle:
        YAML().dump(config, handle)
    record_path = fixture / "fixture.json"
    record = json.loads(record_path.read_text(encoding="utf-8"))
    record.update(
        target_mu_ev=float(target_mu_ev),
        fixed_atom_indices_1based=list(range(1, 9)),
        mobile_atom_indices_1based=list(range(9, 13)),
        md_integrator="fixed-cell Langevin NVT",
    )
    record["changes"].append("Fix bottom eight Pt atoms and leave top four Pt atoms mobile")
    record["sha256"] = {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in fixture.iterdir()
        if path.is_file() and path.name != "fixture.json"
    }
    record_path.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    return path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--target-mu", type=float, required=True)
    parser.add_argument("--steps", type=int, default=5)
    args = parser.parse_args()
    print(materialize_md(args.destination, args.target_mu, steps=args.steps))
