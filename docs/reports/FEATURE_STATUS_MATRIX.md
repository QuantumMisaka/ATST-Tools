# ATST-Tools Feature Status Matrix

**Version**: 2.3.2
**Last Updated**: 2026-10-20
**Status**: 2.3.2 release candidate; publication status is tracked in the release notes
**Owner**: ATST-Tools maintainers

This matrix describes the current source tree's supported, experimental and
unsupported capabilities. Configuration details belong in the
[user reference](../user/CONFIG_REFERENCE.md); source, artifact and
installed-package verification belong in the
[2.3.2 release notes](../releases/RELEASE_NOTES_2.3.2.md).

“Supported” describes implemented workflow behavior. Scientific convergence,
transition-state validation and downstream deployment each require their own
evidence. Constant-potential and DMF retain their experimental scope.

| Feature | Status | Description | Notes |
| :--- | :--- | :--- | :--- |
| **Relax** | ✅ Supported | Geometry Optimization | Uses ASE optimizers. |
| **Constant-potential** | 🧪 Development candidate | ABACUS fixed-geometry single points/serial scans; fixed-cell `relax`/`neb` and ASE `md` with `compensated_gate` (MD added in 2.3.2). | `reference_fcp` is reference-only. Compensated targets use custom `mu` in eV and same-run gate/dipole density facts. No cell-relax/stress, AutoNEB, D2S, native/NPT CP MD, two-Fermi or explicit `nupdown`. [Configuration](../user/CONFIG_REFERENCE.md#211-constant-potential-evaluation-development-candidate); [validation scope](#constant-potential-validation). |
| **Vibration** | ✅ Supported | Frequency Analysis and TS validation | Finite difference method with JSON results, TS validation, and artifact manifest support. |
| **NEB** | ✅ Supported | Nudged Elastic Band | CI-NEB, two-stage NEB, endpoint single-point repair, optional endpoint relaxation, native ASE selector, artifact manifest, ABACUS STRU inputs for `atst neb make`, and MPI image-level parallelism are supported. |
| **AutoNEB** | ✅ Supported | Automated NEB | Adaptive image handling, native ASE selector, endpoint single-point repair, and MPI image-level parallelism are supported. |
| **Dimer** | ✅ Supported | TS Search | Min-mode following. |
| **Sella** | ✅ Supported | Saddle Point Finder | Robust optimization. |
| **CCQN** | ✅ Supported | Constrained Cone Quasi-Newton TS Search | Standalone workflow and D2S refinement option, including reactive-mode enumeration, product alignment, diagnostics, mode manifest, and artifact manifest. |
| **D2S** | ✅ Supported | Double-Ended to Single | `atst run` dispatches endpoint optimization, rough DyNEB, then Dimer, Sella, or CCQN refinement. Experimental `rough_method: dmf` can replace rough DyNEB but is not a supported production default. |
| **IRC** | ✅ Supported | Intrinsic Reaction Coordinate | Sella backend and descent backend are supported, with controlled boundary diagnostics and artifact manifests. |
| **MD** | ✅ Supported | ASE-driven and ABACUS-native molecular dynamics. | Native MD selects the `OUT.<suffix>` directory from INPUT, supports flat outputs, and rejects runs with no frames. Execution and artifact-handoff evidence is recorded in the [integration closeout](ABACUSLITE_230_CLOSEOUT_2026-10-19.md); convergence depends on the calculation. |
| **Independent abacuslite required backend** | ✅ Supported | Installed calculator, input/output readers, reverse configuration, native MD and gate compensation. | Requires `abacuslite>=1.2.1,<2`; missing or broken imports fail. Historical import paths forward to the same installed package. [Wrapper guide](../user/ABACUSLITE_WRAPPER_GUIDE.md); [release verification](../releases/RELEASE_NOTES_2.3.1.md). |
| **ABACUSLite SCF frame identity** | ✅ Supported | Periodic-aware running-log frame selection | The installed `abacuslite==1.2.1` baseline carries the `1e-5 Å`, `rtol=0` SCF identity behavior tested by ATST; the retired vendored copy was the previous baseline. Native `relax`/`md`/`MD_dump` retain last-frame semantics. |
| **Fixed-cell input and NEB runtime preparation** | ✅ Implemented (2.2.9) | Input fidelity and common periodic runtime boundaries | `1e-5 Å`, `rtol=0`; preserve input rounding, warn/override larger cell or PBC differences on copies. NEB/AutoNEB/DyNEB adapts only necessary periodic runtime rows, invalidates changed-boundary endpoint caches, and records actual boundaries in trajectory and manifest. ASE source unchanged. Local controlled-calculator evidence; no new GPU evidence. |
| **DMF** | 🧪 Experimental | Direct MaxFlux TS candidate/path optimizer | Standalone `calculation.type: dmf` and D2S `rough_method: dmf` are available for candidate generation. Outputs are TS candidates, not validated TS results. PBC support is limited to explicit `cartesian_unwrapped` experimental mode. Requires `cyipopt`/IPOPT at runtime. |
| **Artifact Manifests** | ✅ Supported | Workflow output registry | Implemented for NEB, D2S, CCQN, Vibration, IRC, MD, AutoNEB, Sella, Relax, and experimental DMF. |
| **Convergence Stage Records** | ✅ Supported | Optimizer-owned convergence facts in manifest stages | Every optimize workflow records a strict tri-state `converged` (`true`/`false`/`null`) plus stage identity (`role`, `criterion`, IRC `direction`, AutoNEB `subset`, thresholds, stage-local steps). `null` is unknown and never inferred from trajectories or raw force maxima; advisories are diagnostic only and never change completion or exit codes. Documented in `docs/user/PYTHON_API_REFERENCE.md`; contract fixtures in `docs/reports/data/convergence_fixtures_20260917/`. New in 2.2.6; additive stage fields. |
| **English Runtime Diagnostics** | ✅ Supported | Program-visible output language | `print`/warning/logging output, exception messages and CLI help are English across the package, enforced by `tests/unit/test_program_output_language.py`; documentation keeps its existing language and the historical abacuslite compatibility path is excluded. New in 2.2.6; the printed diagnostics language is the only user-visible change. |
| **API Process Runner** | ✅ Supported | External-host API handoff | `python -m atst_tools.api.runner` writes root-only `atst-api-result-v1` JSON and preserves caller-owned scheduler/MPI launch. |
| **NDJSON Progress Events** | ✅ Supported | Structured progress observability | `RunOptions(progress=True)` / runner `--progress` emit one JSON line per event (`workflow_start`, then one `image_step` per NEB/AutoNEB band image) and forward the same mapping to `progress_callback`. |
| **Plotting Helpers** | ✅ Supported | Energy visualization | Stable `neb_energy_profile`, `sella_energy_curve`, `ccqn_energy_curve` API helpers plus the `python -m atst_tools.utils.plot` CLI adapter; matplotlib is an optional `[plot]` extra. |
| **Result Profiles/Plots Extensions** | ✅ Supported | Opt-in result-envelope fields | `RunOptions(profiles=True)` / `--profiles` and `RunOptions(plots=True)` / `--plots` add optional per-image/per-step summaries and plot PNG paths to `atst-api-result-v1` documents without changing the established fields. |
| **Image-Level MPI Parallelism** | ✅ Supported | ASE NEB/AutoNEB image parallelism | Requires MPI-launched Python and compatible `mpi4py`; ABACUS nested MPI remains site-launcher dependent. |
| **Portable MPI Diagnostics** | ✅ Supported | Recovery guidance for missing `mpi4py` in image-parallel launches | Names the explicit `parallel` extra and site-compatible `MPICC` source rebuild; does not promise arbitrary MPI ABI compatibility. New in 2.2.4; additive diagnostic/documentation change. |
| **Runtime Device Binding & Evidence** | ✅ Supported | Device selection, thread budgets and optional run evidence. | `runtime` YAML and CLI/runner controls resolve devices against the caller's binding and trusted allocation. `round_robin` requires a verified pool. Evidence records device/environment facts, counters and sampled memory; unavailable MPS process attribution is explicit. [Runtime configuration](../user/CONFIG_REFERENCE.md#4-runtime-section); [runtime and downstream validation](#runtime-and-downstream-validation). |
| **atst prepare** | ✅ Supported | Reverse config generation from an ABACUS run directory | `atst prepare` and `atst_tools.api.build_config_from_abacus_dir` generate runnable transition YAML (currently NEB) from `INPUT`/`STRU`/`KPT`/`PP`/`ORB`. User-controlled ABACUS values are kept verbatim with three technical floors (`calculation`→`scf`, `cal_force`→`1`, KPT line-mode rejected) and an optional endpoint energy+forces gate (`--no-gate` skips it). New in 2.2.3; backward-compatible. |
| **Trajectory Stress Retention** | ✅ Supported | NEB/AutoNEB trajectories carry per-image stress when the ABACUS calculator requests it | `cal_stress=1` retains per-image stress on NEB/AutoNEB trajectories (serial and image-parallel); `cal_stress=0` behavior is identical to 2.2.x. New in 2.2.3; additive. |
| **GA** | ❌ Not Supported | Genetic Algorithm | ASE 3.28.0 moved GA implementation to the standalone `ase-ga` project; ATST-Tools does not expose GA workflows. |


## Constant-potential validation

The [Pt example](../../examples/19_constant_potential_Pt/README.md) records
bounded energy-derivative and force checks for the compensated boundary,
including numerical and physical-confinement limitations. The
[CP/runtime joint report](ATST_CP_RUNTIME_JOINT_VALIDATION_2026-09-21.md)
records a gate single point and three-point scan with runtime evidence.

The [2026-10-20 parity report](ATST_CP_PARITY_2026-10-20.md) distinguishes
historical LTS `reference_fcp` single-point agreement and develop beta6/beta9/
beta11 execution checks from the corrected release. Those GPU runs used the
pre-review energy mapping. The released abacuslite 1.2.1 preserves ASE free
energy and exposes explicit no-entropy metadata; offline trace replay verifies
the reference convention. Historical GPU runs do not certify the corrected
release or a converged reaction barrier. CP facts-v1 compatibility and the
`profile_status` label are described in the
[configuration reference](../user/CONFIG_REFERENCE.md#31a-constant-potential-decorator).

The 2.3.2 fixed-cell ASE CP MD extension adds per-force electronic
convergence, complete-step checkpoints, validated restart and offline CP
trajectory diagnostics. Its independent analytical tests and dated GPU
observations are recorded in the [MD validation report](ATST_CP_MD_VALIDATION_2026-10-20.md).
The extension retains experimental status; the GPU evidence covers short runs.

## Runtime and downstream validation

P0–P5 and the ATST runtime delivery are complete. The
[SAI V100 report](ATST_RUNTIME_SAI_V100_VALIDATION_2026-09-21.md) records the
DP/ABACUS matrix, multi-card NEB, shared-device checks, 8-image/8-card case and
host/SIF pairing. [P2 closeout](ATST_P2_CLOSEOUT_2026-09-22.md) covers read-only
model sharing, restart attribution and MPS evidence. Site MPI rebinding passed
three repeat runs; see the [active GPU plan](../superpowers/plans/2026-09-20-atst-gpu-node-tuning-plan.md)
for that evidence and the DP container/FT²DP conditional items. Dedicated TF
artifact validation is not planned.

Downstream delivery is tracked in the [P6 takeover report](ATST_GPU_P6_TAKEOVER_2026-10-16.md#9-2026-10-20-下游状态复核).
The parent merged the P6 handoff fixes and deployed Toolbox 1.2.76 to dev/prod,
with GaN relax/band/DOS Agent/Flow acceptance. Same-input 1/4-card NEB,
mixed-parallel exploration and remaining default-path acceptance stay open.
The ATST 2.3.1 / abacuslite 1.2.1 package release does not update that runtime.
Constant-potential Paimon/public tool-chain acceptance is not established.

## Earlier implementation and release evidence

- [2.3.0 integration closeout](ABACUSLITE_230_CLOSEOUT_2026-10-19.md): external
  backend migration, installed consumers and disposition of historical worktrees.
- [2.2.9 release](../releases/RELEASE_NOTES_2.2.9.md): fixed-cell preparation,
  shared CLI/API execution and stored-data TS screening.
- [Architecture/distribution review](ARCHITECTURE_DISTRIBUTION_REVIEW_2026-10-13.md),
  [main-branch review](MAIN_BRANCH_REVIEW_2026-10-15.md) and
  [repair acceptance](MAIN_REVIEW_FIXES_2026-10-16.md): findings, corrections
  and their validation limits.
- [Architecture review](ARCHITECTURE_REVIEW_2026-10-16.md): module ownership
  and runtime/benchmark separation.
- [Release archive](../archive/README.md): earlier versions and dated evidence.
