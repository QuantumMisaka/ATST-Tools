# Constant-potential MD implementation plan

This implementation plan is authorized for execution based on the user's 2026-10-20 request to
add constant-potential MD before further upstream discussion. Execution is complete for the implementation delivery; CP MD is now an
unreleased experimental source capability. Baseline: main `d41eba7`.

## Approach and scope

Extend the existing ASE MD workflow with the existing
`ConstantPotentialCalculator`. At each geometry where forces are required,
converge the electron count to the configured chemical potential before
returning forces. Reuse ABACUS through the installed abacuslite package and
existing file-based evaluation directories. Neither ABACUS PR #8136 nor an
abacuslite socket extension is a prerequisite.

The initial supported combination is:

- `calculation.type: md`, `driver: ase`, fixed cell.
- `calculator.name: abacus`, `energy_boundary: compensated_gate`,
  `reference_electrode: custom`, one scalar `target_mu_ev`.
- Existing fixed-cell integrators: VelocityVerlet under the existing `nve`
  selector; Bussi, Langevin and NVTBerendsen under `nvt`.
- Existing SCF/CP convergence checks, fixed Hamiltonian identity, gate/dipole
  conditions, common-Fermi spin restrictions and nuclear constraints.

Reject native ABACUS MD, NPT/variable cell, target lists and `reference_fcp`
for this workflow. Do not change ordinary MD behavior or add another workflow
name, electronic mass, potentiostat, or relaxed-convergence mode.

The `nve` YAML selector retains its existing meaning as an integrator choice.
With an electron reservoir, this is not a fixed-electron-number NVE ensemble.
Describe the dynamics in terms of the effective grand-potential surface.
NVTBerendsen retains its existing thermostat semantics; its inclusion does not
establish exact canonical sampling.

## Scientific contract

Reuse the compensated boundary:

- `mu_calc = E_Fermi + compensation_derivative`.
- `Omega = E_KohnSham - mu_target * (N - N0)` for the current backend energy
  convention. Preserve force-consistent `free_energy` bookkeeping and do not
  substitute the reference-only no-entropy energy.
- Only a converged SCF and a CP residual within tolerance may supply forces.
- No SCF/CP result from a different geometry or electron count may be reused
  as the current force. A previous electron count is an initial guess only.
- For unthermostatted fixed-target dynamics, diagnose nuclear kinetic energy
  plus the force-consistent grand potential, including any supported nuclear
  constraint contributions once. Raw DFT energy plus kinetic energy is not
  the appropriate conserved quantity when electrons exchange with a reservoir.

The upstream FCP v2 README's MD suggestion (`always_adjust=False`,
`max_FCP_iter=1`) permits potential fluctuations. Its implementation returns
results when the iteration budget expires. ATST currently raises instead, and
this first MD implementation preserves ATST's strict convergence contract.
See the [upstream README](https://github.com/YuLiu98/FCP-abacus-ase/blob/9f96cacb6b69cb39aeb180e73f415f806e16a75e/README.md)
and [source](https://github.com/YuLiu98/FCP-abacus-ase/blob/9f96cacb6b69cb39aeb180e73f415f806e16a75e/version2/FCPelectrochem.py).

## Findings from the current code

- `utils/config_schema.py` and `calculators/factory.py` both reject MD with a
  CP decorator. Both entry points need equivalent scope checks.
- `workflows/md.py:AseMDRunner` already builds calculators through the shared
  factory and delegates integration to ASE. No new integrator is needed.
- `ConstantPotentialCalculator.calculate()` persists validated facts in
  `Atoms.info`, but starts each new evaluation from `initial_electrons`.
- Relax/NEB reuse `constant_potential_restart_initial_electrons_for_config()`;
  MD does not yet use it.
- MD velocity initialization currently treats zero momenta as uninitialized.
  This must not resample a valid zero-momentum CP restart.
- MD passes a trajectory filename without selecting append behavior. ASE
  3.28.0 defaults to overwriting; CP restart must preserve prior trajectory
  frames explicitly.
- Current MD summaries contain energy and temperature but no CP residuals,
  actual completed-step count, or electronic failure record. Current offline
  postprocessing numbers stored frames rather than actual integration steps.

## Task 1: Configuration and factory integration

- [x] Add focused rejection/acceptance tests before changing behavior in
  `tests/unit/test_constant_potential.py` / `tests/unit/test_config.py` and MD
  tests. Cover CLI/API normalization and direct factory/runner use.
- [x] Permit only the combinations above. Validate the same restrictions in
  direct execution paths so native MD cannot silently ignore a CP block.
- [x] Require one scalar target, fixed cell and the existing compensated
  boundary. Keep existing AutoNEB/D2S and other unsupported combinations
  rejected. Reuse validation helpers when checks would otherwise diverge.
- [x] Update schema descriptions and regenerate the YAML variable table.

Acceptance: supported configurations reach the CP factory; invalid combinations
fail before structure mutation, velocity initialization or backend launch.

## Task 2: Integrate complete MD steps and electronic initial guesses

- [x] Add a CP-specific adapter/observer to `AseMDRunner` rather than a second
  workflow or a copy of ASE's integrators. Keep it small; factor cohesive CP
  bookkeeping into a helper if needed.
- [x] Evaluate the initial state and validate its CP facts before recording
  frame zero. Register completed-step handling before ordinary output
  observers, or own the output attachments where ASE ordering requires it.
- [x] After a complete accepted step, carry only that state's validated
  electron count into the next CP solve. Keep the configured neutral reference
  and target unchanged. Do not modify the default electronic iteration trace
  of existing single-point, scan, Relax or NEB workflows. Do not invalidate
  the completed state's cached results just to stage the next initial guess:
  `ConstantPotentialCalculator.set()` currently clears those results.
- [x] Preserve ASE constraints and apply their energy/force contributions once.
  Check that ordinary logging and trajectory writing consume cached results
  without triggering extra CP solves.
- [x] Preserve existing initial-velocity behavior for fresh runs. On CP restart,
  require valid finite momenta and retain them even when all are zero.

Acceptance: a geometry-dependent analytical electronic backend runs through
actual ASE integration and the real CP solver; its converged electron number
and forces agree with the independent analytic solution at every saved step.

## Task 3: Results, restart and failure handling

- [x] Reuse the existing CP facts envelope on saved atoms; add MD step/time
  metadata without duplicating the CP identity or inventing another CP schema.
- [x] Keep `md.traj`, `md.log`, `md_summary.json` and the manifest. Register a
  CP diagnostic JSONL file and a single-frame CP MD checkpoint under the
  configured MD `directory`, with explicit roles in the existing manifest.
  These outputs are automatic for CP MD; no extra output-path YAML fields are
  needed initially.
- [x] Record each completed step's index/time, electron number, target,
  calculated chemical potential, residual, CP iteration count, raw backend
  energy, force-consistent grand potential, nuclear kinetic energy and
  temperature. Label energy definitions in the summary/log metadata.
- [x] Atomically replace the checkpoint after each complete step, including
  positions, cell, constraints, momenta, step/time and validated CP facts.
  Write the same checkpoint for the validated initial state at step zero, so
  failure during the first integration step still leaves a restart point.
  Checkpoints and CP diagnostic rows are written every complete step; regular
  trajectory and MD log output follow the configured interval.
- [x] `restart: true` for CP MD reads that checkpoint and validates geometry,
  assets, target, boundary and residual with the existing restart helper.
  Validate the MD timestep/integrator settings used for continuation. Treat
  `steps` as additional steps; preserve global step/time and append outputs
  without a duplicate initial frame. Anchor regular output cadence to the
  global step, not the resumed segment's local step. With interval 5 and a
  checkpoint at step 7, a regular trajectory/log tail at step 5 is valid, and
  the next scheduled record is step 10, not step 7 or 12. Check output tails
  against their respective cadence and shared step/time identity; reject
  future, mismatched or otherwise inconsistent tails before modifying them.
  Do not silently truncate or mix runs.
- [x] Document restart as continuation of nuclear state and converged electron
  count. It does not restore a partial SCF or guarantee bitwise continuation of
  stochastic thermostats; full random-generator/thermostat checkpoints are
  outside this first implementation.
- [x] On SCF/CP failure, retain the last complete checkpoint and trajectory,
  write a failed summary/manifest with requested and completed steps and the
  available electronic diagnostics, and propagate the failure through CLI/API.
  Never publish the partially advanced nuclear state as a completed frame.
  Close output handles on both success and failure.
- [x] Extend `utils/md_post.py` to read CP data from saved frames without a
  calculator. Validate persisted facts before displaying them. Use saved MD
  step/time when available and retain the ordinary-MD fallback for old files.

Acceptance: interrupted and resumed deterministic trajectories match a
continuous reference within numerical tolerance; zero momenta survive restart;
a split at step 7 with interval 5 preserves checkpoint progress and output
cadence; first-step failure preserves the validated step-zero checkpoint. Failed CP
steps do not appear as successful MD frames or manifests.

## Task 4: Numerical and regression verification

- [x] Use an analytic backend with geometry/electron coupling, for example
  `E(x,n)=k*x^2/2+(n-n0-a*x)^2/(2*C)` with positive `k,C`.
  Its independent stationary solution is `n*=n0+a*x+C*mu_target`, and the
  relaxed grand-potential force is `-k*x+a*mu_target`. Exercise the production
  CP loop with synthetic, internally consistent compensation facts.
- [x] Test force versus finite differences of the minimized grand potential,
  successful and deliberately failed electronic solves, metadata freshness,
  constrained atoms, extra logging reads, and one-atom/free-DOF edge cases
  relevant to the chosen thermostat.
- [x] Test fixed-cell NVT integrators and use VelocityVerlet to compare energy
  drift and timestep refinement. Set tolerances from analytic error behavior,
  rather than using an arbitrary universal DFT acceptance threshold.
- [x] Test restart with missing/stale CP facts, changed target/assets/cell,
  absent/nonfinite momenta, output-tail conflicts and failure during a step.
- [x] Run the existing CP/reference/restart, gate, MD/postprocessing, config,
  CLI/API and example suites. Run the full offline suite after integration.

## Task 5: Real ABACUS example and evidence

- [x] Extend `examples/19_constant_potential_Pt/` with a materialized CP MD
  configuration and instructions. Reuse its PP/ORB and explicit gate boundary;
  use a fixture with enough mobile degrees of freedom for the chosen thermostat.
  Do not transfer a target chemical potential from a different Hamiltonian.
- [x] Use Slurm on the project-required SAI 4V100 resources with LCAO
  `ks_solver: cusolver`, an isolated run directory and the installed external
  abacuslite baseline. Record the exact engine/package environment.
- [x] First compare force and grand-potential finite differences for the chosen
  moving geometry, then run a short fixed-cell NVT trajectory, a split/restart
  run and an unthermostatted timestep comparison. Record per-step SCF/CP
  convergence, electron counts, forces, energies and wall time. Keep electronic
  tolerance refinement separate from integration timestep refinement.
- [x] Inspect gate/block geometry and charge confinement throughout the run.
  Existing single-point or force checks do not establish long-time MD validity.
  An unavailable site run remains pending evidence, not a simulated success.

## Task 6: Documentation, review and integration

- [x] Update `CONFIG_REFERENCE.md`, generated YAML documentation, the Chinese
  user guide, examples index, MD/API result descriptions where affected, and
  `FEATURE_STATUS_MATRIX.md` with exactly the tested scope.
- [x] Add a dated evidence report for real runs. Record the feature under the
  next release's notes when its version is selected; do not retroactively
  claim the published 2.3.1 artifacts include it.
- [x] Run documentation governance and relevant checks; obtain independent
  final review of scientific semantics, restart/output integrity and backward
  compatibility before integration. Merge only the checked revision.
- [x] After delivery, absorb the lasting contracts into the user/developer
  references and archive this completed plan under the existing docs rules.

Socket integration remains a separate improvement to the backend transport.
This implementation should provide a useful file-based numerical baseline for
that later change, rather than depending on unmerged engine code.


## Execution record — 2026-10-20

- Workflow: L3, subagent-driven development, isolated worktree
  `feat/constant-potential-md-20261020`, BASE
  `d41eba7a952bc0de381838a87fd4efc6e088100e`.
- Before dispatch all tracked paths were clean; the only untracked file was
  this plan. Package changes are attributable against BASE; no user code edits
  were present in the implementation worktree.
- Baseline: 19 existing MD, MD postprocessing and CP restart tests passed with
  Python 3.10 / ASE 3.28.0 in the retained verification environment.
- Configuration package: 118 existing/new configuration and CP tests passed;
  generated YAML reference and documentation governance passed.
- Runtime/restart, analytic verification and documentation packages completed.
  Full offline suite: 1450 passed, 37 skipped, 1 warning. Two independent
  code-review findings were fixed and re-reviewed successfully.
- Real SAI validation: jobs 1872958 and 1873489 completed with exit 0:0 on
  4V100/LCAO/cusolver using abacuslite 1.2.1. The latter used the final reviewed
  source and verified short MD plus split restart. Force differences, timestep
  comparisons, charge confinement and numerical limitations are retained in
  the dated MD validation report.
- Ruling: deliver the reviewed implementation on the isolated feature branch;
  merge, release-version selection, release notes and publication remain a
  later integration action. This execution request does not require them.
  The conditional release-note/merge steps above are therefore satisfied by
  retaining this explicit boundary, not by claiming an unpublished release.
- Ruling: real-run timestep halving and electronic tolerance variation are
  recorded observations. The short 0.5 fs trajectory does not establish a
  long-time sampling or production-accuracy claim.
- Release version remains 2.3.1; this branch is unreleased implementation work.
  Published artifacts are not being relabeled or replaced.

Implementation evidence: [dated MD validation report](../../reports/ATST_CP_MD_VALIDATION_2026-10-20.md).
