"""Run one FCP v2 reference calculation and write a normalized record.

This maintainer tool deliberately keeps the FCP v2 source external: the
reference implementation is GPL-3.0 and is not vendored into ATST-Tools.
Use it only in an environment that provides the legacy ``ase-abacus`` package
and an ABACUS executable.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import traceback
from typing import Any, Mapping

import numpy as np
from ase.io import read


_REQUIRED_TOP_LEVEL = (
    "atoms",
    "command",
    "pseudo_dir",
    "pseudopotentials",
    "basissets",
    "parameters",
    "fcp",
)
_REQUIRED_FCP = ("U", "NELECT", "NELECT0")


def _require_keys(mapping: Mapping[str, Any], keys: tuple[str, ...], context: str) -> None:
    """Raise a clear error when a required configuration key is absent."""
    missing = [key for key in keys if key not in mapping]
    if missing:
        raise ValueError(f"{context} is missing required keys: {', '.join(missing)}")


def _load_fcp_class(source: Path):
    """Load the external FCP v2 source and return its ``FCP`` class."""
    spec = importlib.util.spec_from_file_location("fcp_v2_reference", source)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load FCP v2 source: {source}")
    module = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(module)
    except ImportError as exc:
        raise RuntimeError(
            "FCP v2 reference requires the legacy ase-abacus package and pandas; "
            f"import failed: {exc}"
        ) from exc
    if not hasattr(module, "FCP"):
        raise RuntimeError(f"FCP v2 source does not define FCP: {source}")
    return module.FCP


def _read_last_fcp_log(path: Path) -> dict[str, float | int | None]:
    """Read the last evaluated state from FCP v2's ``tmp-log-FCP.txt``."""
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines()]
    lines = [line for line in lines if line]
    if len(lines) < 2:
        raise ValueError(f"FCP log has no data rows: {path}")
    header = lines[0].split("\t")
    index = {name: position for position, name in enumerate(header)}
    required = (
        "loop",
        "NELECT",
        "Fermi(eV)",
        "Fermishift(eV)",
        "mu(eV)",
        "conv(V)",
        "Etoten(eV)",
        "Etoten_grand(eV)",
        "Cpersurf(e/V/A^2)",
    )
    missing = [name for name in required if name not in index]
    if missing:
        raise ValueError(f"FCP log is missing columns: {', '.join(missing)}")
    rows = [line.split("\t") for line in lines[1:]]
    row = rows[-1]
    if len(row) <= max(index[name] for name in required):
        raise ValueError(f"FCP log last row is incomplete: {path}")

    def value(name: str) -> float | None:
        try:
            return float(row[index[name]])
        except (TypeError, ValueError):
            return None

    return {
        "iterations": int(float(row[index["loop"]])),
        "nelec": value("NELECT"),
        "efermi_eV": value("Fermi(eV)"),
        "fermishift_eV": value("Fermishift(eV)"),
        "mu_calc_eV": value("mu(eV)"),
        "residual_mu_eV": value("conv(V)"),
        "raw_energy_eV": value("Etoten(eV)"),
        "energy_eV": value("Etoten_grand(eV)"),
        "capacitance_per_area": value("Cpersurf(e/V/A^2)"),
    }


def _build_record(
    *,
    log: Mapping[str, Any],
    fcp: Any,
    atoms: Any,
    config: Mapping[str, Any],
    source: Path,
    source_sha256: str,
    energy: float,
    free_energy: float,
    forces: np.ndarray,
) -> dict[str, Any]:
    """Build the comparator record from one completed FCP v2 evaluation."""
    fcp_parameters = dict(config["fcp"])
    mu_target = -(
        float(fcp_parameters["U"]) + float(fcp_parameters.get("work_ref", 4.6))
    )
    fermishift = log.get("fermishift_eV")
    vacuum = None if fermishift is None else -float(fermishift)
    capacitance = getattr(fcp, "C", None)
    if capacitance is None and log.get("capacitance_per_area") is not None:
        area = float(np.linalg.norm(np.cross(atoms.cell[0], atoms.cell[1])))
        capacitance = float(log["capacitance_per_area"]) * area
    residual = log.get("residual_mu_eV")
    tolerance = float(fcp_parameters.get("FCPconv", 0.01))
    converged = None if residual is None else abs(float(residual)) <= tolerance

    return {
        "nelec": log.get("nelec"),
        "efermi_eV": log.get("efermi_eV"),
        "vacuum_level_eV": vacuum,
        "fermishift_eV": fermishift,
        "mu_calc_eV": log.get("mu_calc_eV"),
        "mu_target_eV": mu_target,
        "residual_mu_eV": residual,
        "energy_eV": float(energy),
        "free_energy_eV": float(free_energy),
        "raw_energy_eV": log.get("raw_energy_eV"),
        "forces_eV_per_A": np.asarray(forces, dtype=float).tolist(),
        "capacitance_e_per_V": None if capacitance is None else float(capacitance),
        "iterations": log.get("iterations"),
        "converged": converged,
        "reference_source": str(source),
        "reference_sha256": source_sha256,
        "fcp_parameters": fcp_parameters,
        "engine_command": config["command"],
        "directory": str(config.get("directory", "fcp_reference_run")),
    }


def run_reference(config_path: Path, source_path: Path) -> dict[str, Any]:
    """Run one FCP v2 calculation from a JSON configuration.

    Args:
        config_path: JSON configuration describing atoms, ABACUS, and FCP.
        source_path: Path to the pinned ``version2/FCPelectrochem.py``.

    Returns:
        The normalized reference record.

    Raises:
        RuntimeError: If ase-abacus, pandas, or the FCP source cannot be loaded.
        ValueError: If the configuration is incomplete.
    """
    config = json.loads(config_path.read_text(encoding="utf-8"))
    _require_keys(config, _REQUIRED_TOP_LEVEL, "FCP reference configuration")
    _require_keys(config["fcp"], _REQUIRED_FCP, "FCP reference fcp block")
    if "basis_dir" not in config and "orbital_dir" not in config:
        raise ValueError(
            "FCP reference configuration needs basis_dir or orbital_dir"
        )

    try:
        from ase.calculators.abacus import Abacus, AbacusProfile
    except ImportError as exc:
        raise RuntimeError(
            "FCP v2 reference requires the legacy ase-abacus package"
        ) from exc

    source_path = source_path.expanduser().resolve()
    source_sha256 = hashlib.sha256(source_path.read_bytes()).hexdigest()
    fcp_class = _load_fcp_class(source_path)

    directory = Path(config.get("directory", "fcp_reference_run")).expanduser().resolve()
    directory.mkdir(parents=True, exist_ok=True)
    atoms_path = Path(config["atoms"])
    atoms_format = config.get("atoms_format")
    if atoms_format is None and atoms_path.name.upper() == "STRU":
        atoms_format = "abacus"
    atoms = read(atoms_path, format=atoms_format)

    profile = AbacusProfile(command=config["command"])
    inner_parameters = dict(config["parameters"])
    inner_parameters.update(
        {
            "pp": config["pseudopotentials"],
            "basis": config["basissets"],
            "pseudo_dir": config["pseudo_dir"],
            "basis_dir": config.get("basis_dir", config.get("orbital_dir")),
        }
    )
    if "kpts" in config:
        inner_parameters["kpts"] = config["kpts"]
    inner = Abacus(
        profile=profile,
        directory=str(directory),
        **inner_parameters,
    )

    fcp_parameters = dict(config["fcp"])
    fcp_parameters.setdefault("C", 1.0 / 80.0)
    fcp_parameters.setdefault("FCPconv", 0.01)
    fcp_parameters.setdefault("work_ref", 4.6)
    fcp_parameters.setdefault("max_FCP_iter", 100)
    fcp_parameters.setdefault("always_adjust", True)
    fcp_parameters.setdefault("adaptive_lr", False)
    fcp_parameters.setdefault("FCPmethod", "Newton-fitting")
    fcptxt = Path(fcp_parameters.pop("fcptxt", directory / "log-fcp.txt")).expanduser()
    fcptxt.parent.mkdir(parents=True, exist_ok=True)

    fcp = fcp_class(innercalc=inner, fcptxt=str(fcptxt), **fcp_parameters)
    fcp.directory = str(directory)
    atoms.calc = fcp
    energy = float(atoms.get_potential_energy())
    free_energy = float(atoms.get_potential_energy(force_consistent=True))
    forces = atoms.get_forces()

    log = _read_last_fcp_log(directory / "tmp-log-FCP.txt")
    return _build_record(
        log=log,
        fcp=fcp,
        atoms=atoms,
        config={**config, "directory": str(directory)},
        source=source_path,
        source_sha256=source_sha256,
        energy=energy,
        free_energy=free_energy,
        forces=forces,
    )


def main(argv: list[str] | None = None) -> int:
    """Run the command-line reference producer."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="FCP reference JSON config.")
    parser.add_argument(
        "--fcp-source",
        required=True,
        help="Path to the pinned version2/FCPelectrochem.py.",
    )
    parser.add_argument("--output", required=True, help="Normalized record JSON.")
    args = parser.parse_args(argv)

    try:
        record = run_reference(Path(args.config), Path(args.fcp_source))
    except Exception as exc:
        print(
            json.dumps(
                {
                    "ok": False,
                    "error": str(exc),
                    "traceback": traceback.format_exc(),
                },
                indent=2,
            )
        )
        return 2
    Path(args.output).write_text(
        json.dumps(record, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(record, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
