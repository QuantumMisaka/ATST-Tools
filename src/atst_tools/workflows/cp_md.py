"""Strict constant-potential bookkeeping for the ASE molecular-dynamics runner."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from pathlib import Path
from typing import Any, Callable

import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.io.trajectory import Trajectory

from atst_tools.calculators.constant_potential import (
    ConstantPotentialCalculator,
    constant_potential_restart_initial_electrons_for_config,
    read_constant_potential_facts,
)
from atst_tools.calculators.factory import CalculatorFactory
from atst_tools.utils.artifacts import write_artifact_manifest
from atst_tools.utils.md_post import run_md_workflow_postprocess


STEP_KEY = "atst_md_step"
TIME_KEY = "atst_md_time_fs"
SETTINGS_KEY = "atst_md_settings"
LOG_HEADER = (
    "# step time_fs grand_potential_eV raw_energy_eV kinetic_energy_eV "
    "temperature_K nelec mu_target_ev mu_calc_ev residual_mu_ev cp_iterations\n"
)


def _settings(config: dict[str, Any]) -> dict[str, Any]:
    """Return the MD controls whose change would alter continuation or cadence."""
    keys = (
        "ensemble", "algorithm", "timestep_fs", "temperature_K", "taut_fs",
        "friction_fs_inv", "force_temperature", "stationary", "zero_rotation",
        "loginterval",
    )
    return {key: config.get(key) for key in keys}


def _step_time(atoms: Atoms, timestep_fs: float) -> tuple[int, float]:
    """Read and validate durable global MD coordinates."""
    step = atoms.info.get(STEP_KEY)
    time_fs = atoms.info.get(TIME_KEY)
    if isinstance(step, bool) or not isinstance(step, (int, np.integer)) or step < 0:
        raise ValueError("CP MD output has an invalid global step")
    try:
        time_value = float(time_fs)
    except (TypeError, ValueError) as exc:
        raise ValueError("CP MD output has an invalid global time") from exc
    if not np.isfinite(time_value) or not np.isclose(
        time_value, int(step) * timestep_fs, rtol=1e-10, atol=1e-10
    ):
        raise ValueError("CP MD output has an inconsistent global step/time")
    return int(step), time_value


def _validate_tail(path: Path, expected_step: int, timestep_fs: float) -> Atoms:
    """Check the last regular frame before opening a trajectory for append."""
    if not path.is_file():
        raise ValueError(f"CP MD restart requires existing trajectory: {path}")
    try:
        tail = read(path, index=-1)
    except Exception as exc:
        raise ValueError("CP MD trajectory tail is unreadable") from exc
    step, _ = _step_time(tail, timestep_fs)
    if step != expected_step or read_constant_potential_facts(tail) is None:
        raise ValueError("CP MD trajectory tail conflicts with checkpoint or output cadence")
    return tail


def _validate_log_tail(path: Path, expected_step: int, timestep_fs: float) -> dict[str, float]:
    """Check the ordinary MD log tail against the checkpoint cadence."""
    if not path.is_file():
        raise ValueError(f"CP MD restart requires existing log: {path}")
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or lines[0] != LOG_HEADER.rstrip("\n") or len(lines) < 2:
        raise ValueError("CP MD log header or tail is invalid")
    fields = lines[-1].split()
    if len(fields) != 11:
        raise ValueError("CP MD log tail is malformed")
    try:
        step = int(fields[0])
        time_fs = float(fields[1])
    except (IndexError, TypeError, ValueError) as exc:
        raise ValueError("CP MD log tail is malformed") from exc
    if step != expected_step or not np.isclose(time_fs, step * timestep_fs, atol=1e-10):
        raise ValueError("CP MD log tail conflicts with checkpoint or output cadence")
    names = (
        "grand_potential_eV", "raw_energy_eV", "kinetic_energy_eV",
        "temperature_K", "nelec", "mu_target_ev", "mu_calc_ev",
        "residual_mu_ev", "cp_iterations",
    )
    try:
        values = {name: float(value) for name, value in zip(names, fields[2:])}
    except ValueError as exc:
        raise ValueError("CP MD log tail has invalid numeric fields") from exc
    if not all(np.isfinite(value) for value in values.values()):
        raise ValueError("CP MD log tail has non-finite numeric fields")
    return values


def _validate_diagnostics_tail(path: Path, checkpoint_step: int, timestep_fs: float) -> dict[str, Any]:
    """Require one durable diagnostic row for every completed global step."""
    if not path.is_file():
        raise ValueError(f"CP MD restart requires existing diagnostics: {path}")
    try:
        rows = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError("CP MD diagnostic tail is malformed") from exc
    if len(rows) != checkpoint_step + 1:
        raise ValueError("CP MD diagnostic tail conflicts with checkpoint")
    for step, row in enumerate(rows):
        if row.get("step") != step:
            raise ValueError("CP MD diagnostics have a missing or duplicated step")
        try:
            time_fs = float(row["time_fs"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError("CP MD diagnostic time is malformed") from exc
        if not np.isfinite(time_fs) or not np.isclose(time_fs, step * timestep_fs, atol=1e-10):
            raise ValueError("CP MD diagnostic step/time is inconsistent")
    return rows[-1]


def _same_diagnostic_value(left: Any, right: Any) -> bool:
    """Compare rounded log values against full-precision durable records."""
    try:
        a, b = float(left), float(right)
    except (TypeError, ValueError):
        return False
    return bool(np.isfinite(a) and np.isfinite(b) and np.isclose(a, b, rtol=1e-10, atol=1e-10))


def _nuclear_state_fingerprint(atoms: Atoms) -> str:
    """Fingerprint the saved nuclear state and its ASE constraints."""
    constraints = []
    for constraint in atoms.constraints:
        encoder = getattr(constraint, "todict", None)
        constraints.append({
            "class": f"{type(constraint).__module__}.{type(constraint).__qualname__}",
            "data": encoder() if callable(encoder) else repr(constraint),
        })
    state = {
        "numbers": atoms.numbers.tolist(),
        "positions": atoms.positions.tolist(),
        "cell": atoms.cell.array.tolist(),
        "pbc": atoms.pbc.tolist(),
        "masses": atoms.get_masses().tolist(),
        "momenta": atoms.get_momenta().tolist(),
        "constraints": constraints,
    }

    def _encode(value: Any) -> Any:
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, np.generic):
            return value.item()
        raise TypeError(f"CP MD nuclear state contains unsupported type: {type(value).__name__}")

    serialized = json.dumps(state, sort_keys=True, separators=(",", ":"), default=_encode)
    return hashlib.sha256(serialized.encode("utf-8")).hexdigest()


def _validate_saved_nuclear_row(atoms: Atoms, row: dict[str, Any], source: str) -> None:
    """Bind a diagnostic row to its saved ASE frame without backend calls."""
    if row.get("nuclear_state_fingerprint") != _nuclear_state_fingerprint(atoms):
        raise ValueError(f"CP MD {source} nuclear state disagrees with diagnostic row")
    try:
        grand_potential = float(atoms.get_potential_energy(force_consistent=True))
    except Exception as exc:
        raise ValueError(f"CP MD {source} grand potential cannot be read") from exc
    kinetic = float(atoms.get_kinetic_energy())
    for key, value in (
        ("grand_potential_eV", grand_potential),
        ("kinetic_energy_eV", kinetic),
        ("conserved_quantity_eV", grand_potential + kinetic),
    ):
        if not _same_diagnostic_value(row.get(key), value):
            label = key.replace("_eV", "").replace("_", " ")
            raise ValueError(f"CP MD {source} {label} disagrees with diagnostic row")


def _atomic_checkpoint(path: Path, atoms: Atoms) -> None:
    """Replace a single-frame ASE checkpoint only after a complete step."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, suffix=".traj", delete=False) as handle:
        temporary = Path(handle.name)
    try:
        write(temporary, atoms, format="traj")
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _row(atoms: Atoms, step: int, timestep_fs: float) -> dict[str, Any]:
    """Measure one complete CP state using cached calculator results."""
    payload = read_constant_potential_facts(atoms)
    if payload is None:
        raise ValueError("CP MD step has no valid converged electronic facts")
    facts = payload["facts"]
    kinetic = float(atoms.get_kinetic_energy())
    # ASE applies constraint energy once here; the CP facts record the raw
    # unconstrained electronic contribution from the same cached solve.
    grand_potential = float(atoms.get_potential_energy(force_consistent=True))
    return {
        "step": step,
        "time_fs": step * timestep_fs,
        "nelec": float(facts["nelec"]),
        "mu_target_ev": float(facts["mu_target"]),
        "mu_calc_ev": float(facts["mu_calc"]),
        "residual_mu_ev": float(facts["residual_mu"]),
        "cp_iterations": int(facts["cp_iterations"]),
        "raw_energy_eV": float(facts["raw_energy"]),
        "grand_potential_eV": grand_potential,
        "kinetic_energy_eV": kinetic,
        "temperature_K": float(atoms.get_temperature()),
        "conserved_quantity_eV": grand_potential + kinetic,
        "nuclear_state_fingerprint": _nuclear_state_fingerprint(atoms),
    }


class ConstantPotentialMD:
    """Own CP MD persistence around the existing ASE integrator steps."""

    def __init__(
        self,
        config: dict[str, Any],
        calc_config: dict[str, Any],
        build_dynamics: Callable[[Atoms], Any],
        initialize_velocities: Callable[[Atoms], None],
    ) -> None:
        self.config = config
        self.calc_config = calc_config
        self.build_dynamics = build_dynamics
        self.initialize_velocities = initialize_velocities
        self.directory = Path(calc_config["directory"])
        self.checkpoint = self.directory / "cp_md_checkpoint.traj"
        self.diagnostics_path = self.directory / "cp_md_diagnostics.jsonl"
        self.trajectory_path = Path(calc_config["trajectory"])
        self.log_path = Path(calc_config["logfile"])
        self.interval = int(calc_config["loginterval"])
        self.timestep_fs = float(calc_config["timestep_fs"])

    def _read_restart(self) -> tuple[Atoms, int, float, dict[str, Any]]:
        """Validate the checkpoint and all output tails before any append."""
        if not self.checkpoint.is_file():
            raise ValueError(f"CP MD restart checkpoint is missing: {self.checkpoint}")
        atoms = read(self.checkpoint)
        if atoms.info.get(SETTINGS_KEY) != _settings(self.calc_config):
            raise ValueError("CP MD restart integrator/timestep/output settings do not match checkpoint")
        step, time_fs = _step_time(atoms, self.timestep_fs)
        if not atoms.has("momenta") or not np.isfinite(atoms.get_momenta()).all():
            raise ValueError("CP MD restart requires finite saved momenta")
        nelec = constant_potential_restart_initial_electrons_for_config(atoms, self.config)
        if nelec is None:
            raise ValueError("CP MD restart has no validated electron count")
        expected_tail = (step // self.interval) * self.interval
        trajectory_tail = _validate_tail(self.trajectory_path, expected_tail, self.timestep_fs)
        log_tail = _validate_log_tail(self.log_path, expected_tail, self.timestep_fs)
        last_row = _validate_diagnostics_tail(self.diagnostics_path, step, self.timestep_fs)
        _validate_saved_nuclear_row(atoms, last_row, "checkpoint")
        facts = read_constant_potential_facts(atoms)
        assert facts is not None
        if not _same_diagnostic_value(last_row.get("nelec"), nelec) or not _same_diagnostic_value(
            last_row.get("residual_mu_ev"), facts["facts"]["residual_mu"]
        ):
            raise ValueError("CP MD diagnostic tail disagrees with checkpoint electronic facts")
        if not _same_diagnostic_value(last_row.get("kinetic_energy_eV"), atoms.get_kinetic_energy()) or not _same_diagnostic_value(
            last_row.get("temperature_K"), atoms.get_temperature()
        ):
            raise ValueError("CP MD diagnostic tail disagrees with checkpoint nuclear state")
        tail_facts = read_constant_potential_facts(trajectory_tail)
        assert tail_facts is not None
        if tail_facts["identity"] != facts["identity"]:
            raise ValueError("CP MD trajectory tail belongs to a different electronic calculation")
        if expected_tail == step and tail_facts["geometry_fingerprint"] != facts["geometry_fingerprint"]:
            raise ValueError("CP MD trajectory tail geometry disagrees with checkpoint")
        rows = self.diagnostics_path.read_text(encoding="utf-8").splitlines()
        tail_row = json.loads(rows[expected_tail])
        _validate_saved_nuclear_row(trajectory_tail, tail_row, "trajectory tail")
        for key, expected in (
            ("nelec", tail_facts["facts"]["nelec"]),
            ("residual_mu_ev", tail_facts["facts"]["residual_mu"]),
            ("mu_calc_ev", tail_facts["facts"]["mu_calc"]),
            ("mu_target_ev", tail_facts["facts"]["mu_target"]),
            ("raw_energy_eV", tail_facts["facts"]["raw_energy"]),
            ("kinetic_energy_eV", trajectory_tail.get_kinetic_energy()),
            ("temperature_K", trajectory_tail.get_temperature()),
        ):
            if not _same_diagnostic_value(tail_row.get(key), expected):
                raise ValueError("CP MD trajectory tail disagrees with diagnostic row")
        for key, value in log_tail.items():
            if not _same_diagnostic_value(tail_row.get(key), value):
                raise ValueError("CP MD log tail disagrees with diagnostic row")
        for key in ("nelec", "residual_mu_ev", "mu_calc_ev", "mu_target_ev", "raw_energy_eV"):
            expected = {
                "nelec": "nelec", "residual_mu_ev": "residual_mu", "mu_calc_ev": "mu_calc",
                "mu_target_ev": "mu_target", "raw_energy_eV": "raw_energy",
            }[key]
            if not _same_diagnostic_value(last_row.get(key), facts["facts"][expected]):
                raise ValueError("CP MD diagnostic tail disagrees with checkpoint electronic facts")
        if expected_tail == step and not np.allclose(
            trajectory_tail.get_momenta(), atoms.get_momenta(), rtol=1e-10, atol=1e-12
        ):
            raise ValueError("CP MD trajectory tail momenta disagree with checkpoint")
        return atoms, step, nelec, last_row

    def _validate_fresh_outputs(self) -> None:
        """Refuse to overwrite a prior CP MD run implicitly."""
        paths = (
            self.checkpoint, self.diagnostics_path, self.trajectory_path,
            self.log_path, Path(self.calc_config["summary_file"]),
            Path(self.calc_config["final_structure"]),
            Path(self.calc_config["artifact_manifest"]),
        )
        existing = [str(path) for path in paths if path.exists()]
        if existing:
            raise ValueError("CP MD fresh run would overwrite existing outputs: " + ", ".join(existing))

    def _record(self, atoms: Atoms, row: dict[str, Any], *, regular: bool, fresh: bool) -> None:
        """Persist a complete step before exposing it as a restart point."""
        step = row["step"]
        atoms.info[STEP_KEY] = step
        atoms.info[TIME_KEY] = row["time_fs"]
        atoms.info[SETTINGS_KEY] = _settings(self.calc_config)
        with self.diagnostics_path.open("w" if fresh else "a", encoding="utf-8") as handle:
            handle.write(json.dumps(row, allow_nan=False) + "\n")
            handle.flush()
        _atomic_checkpoint(self.checkpoint, atoms)
        if regular:
            with Trajectory(self.trajectory_path, mode="w" if fresh else "a", atoms=atoms) as trajectory:
                trajectory.write(atoms)
            with self.log_path.open("w" if fresh else "a", encoding="utf-8") as handle:
                if fresh:
                    handle.write(LOG_HEADER)
                handle.write(
                    f"{step} {row['time_fs']:.12g} {row['grand_potential_eV']:.12g} "
                    f"{row['raw_energy_eV']:.12g} {row['kinetic_energy_eV']:.12g} "
                    f"{row['temperature_K']:.12g} {row['nelec']:.12g} "
                    f"{row['mu_target_ev']:.12g} {row['mu_calc_ev']:.12g} "
                    f"{row['residual_mu_ev']:.12g} {row['cp_iterations']}\n"
                )

    def _summary(self, row: dict[str, Any] | None, *, start: int, status: str,
                 error: Exception | None = None, calculator: ConstantPotentialCalculator | None = None) -> dict[str, Any]:
        """Describe requested and completed work without reading partial atoms."""
        completed = row["step"] if row is not None else start
        summary: dict[str, Any] = {
            "workflow": "md",
            "driver": "ase",
            "ensemble": self.calc_config["ensemble"],
            "algorithm": self.calc_config["algorithm"],
            "status": status,
            "steps": self.calc_config["steps"],
            "requested_steps": self.calc_config["steps"],
            "completed_steps": completed,
            "segment_completed_steps": completed - start,
            "final_step": completed,
            "final_time_fs": None if row is None else row["time_fs"],
            "final_energy_eV": None if row is None else row["grand_potential_eV"],
            "final_temperature_K": None if row is None else row["temperature_K"],
            "energy_definition": "force_consistent_grand_potential_plus_nuclear_constraint_energy",
            "conserved_quantity_definition": "grand_potential_plus_nuclear_kinetic_energy",
            "trajectory": str(self.trajectory_path),
            "final_structure": self.calc_config["final_structure"],
            "checkpoint": str(self.checkpoint),
            "diagnostics": str(self.diagnostics_path),
        }
        if error is not None:
            summary["error_type"] = type(error).__name__
            summary["error"] = str(error)
            if calculator is not None:
                summary["electronic_diagnostics"] = calculator.diagnostics
        return summary

    def _finish(self, summary: dict[str, Any], *, postprocess_artifacts: list[dict[str, Any]] | None = None) -> None:
        """Publish a summary and manifest for complete or failed segments."""
        summary_path = Path(self.calc_config["summary_file"])
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_path.write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
        artifacts = [
            {"role": "trajectory", "path": str(self.trajectory_path)},
            {"role": "log", "path": str(self.log_path)},
            {"role": "summary", "path": str(summary_path)},
            {"role": "cp_md_checkpoint", "path": str(self.checkpoint)},
            {"role": "cp_md_diagnostics", "path": str(self.diagnostics_path)},
        ]
        if summary["status"] == "complete":
            artifacts.append({"role": "final_structure", "path": self.calc_config["final_structure"]})
        artifacts.extend(postprocess_artifacts or [])
        write_artifact_manifest(
            self.calc_config["artifact_manifest"],
            workflow="md",
            artifacts=artifacts,
            stages=[{
                "name": "ase_cp_md", "status": summary["status"],
                "requested_steps": summary["requested_steps"],
                "completed_steps": summary["completed_steps"],
            }],
            metadata=summary,
        )

    def run(self, read_initial: Callable[[], Atoms], calc_name: str) -> None:
        """Integrate strict CP states and preserve global restart coordinates."""
        restart = bool(self.calc_config.get("restart"))
        if restart:
            atoms, start, warmstart, last_row = self._read_restart()
        else:
            self._validate_fresh_outputs()
            atoms = read_initial()
            start, warmstart = 0, None
            last_row = None
        self.directory.mkdir(parents=True, exist_ok=True)
        self.trajectory_path.parent.mkdir(parents=True, exist_ok=True)
        self.log_path.parent.mkdir(parents=True, exist_ok=True)
        calculator = CalculatorFactory.get_calculator(
            calc_name,
            self.config,
            directory=self.calc_config["directory"],
        )
        if not isinstance(calculator, ConstantPotentialCalculator):
            raise TypeError("CP MD requires a ConstantPotentialCalculator")
        atoms.calc = calculator
        if restart:
            calculator.stage_initial_electrons(warmstart)
        else:
            # Fresh runs retain the ordinary MD velocity initialization policy.
            self.initialize_velocities(atoms)
        dynamics = self.build_dynamics(atoms)
        dynamics.nsteps = start
        try:
            if not restart:
                atoms.get_forces(md=True)
                last_row = _row(atoms, 0, self.timestep_fs)
                self._record(atoms, last_row, regular=True, fresh=True)
                calculator.stage_initial_electrons(last_row["nelec"])
            else:
                # The saved row retains constraint contributions if the first
                # resumed integration attempt fails before a new complete step.
                assert last_row is not None
            for _ in range(int(self.calc_config["steps"])):
                dynamics.step()
                atoms.get_forces(md=True)
                step = dynamics.nsteps + 1
                row = _row(atoms, step, self.timestep_fs)
                dynamics.nsteps = step
                self._record(atoms, row, regular=step % self.interval == 0, fresh=False)
                last_row = row
                calculator.stage_initial_electrons(row["nelec"])
            write(self.calc_config["final_structure"], atoms)
            summary = self._summary(last_row, start=start, status="complete")
            postprocess = run_md_workflow_postprocess(
                self.calc_config,
                metadata={"driver": "ase", "algorithm": self.calc_config["algorithm"]},
            )
            self._finish(summary, postprocess_artifacts=postprocess["artifacts"])
        except Exception as exc:
            self._finish(self._summary(last_row, start=start, status="failed", error=exc,
                                       calculator=calculator))
            raise
