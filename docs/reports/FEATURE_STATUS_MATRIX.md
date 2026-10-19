# ATST-Tools Feature Status Matrix

**Version**: 2.3.1
**Last Updated**: 2026-10-20
**Status**: Source candidate (unpublished)
**Owner**: ATST-Tools maintainers

Release 2.3.0 made the independently installed `abacuslite` package the sole
ABACUS ASE backend. 2.3.1 requires `abacuslite>=1.2.1` for the entropy-aware
`free_energy` handoff. The vendored snapshot is retired as a runtime fallback;
historical imports forward to the installed package. The latest published
package remains 2.3.0.

Current distribution acceptance and development-branch disposition are tracked
in the [2.3.0 integration closeout](ABACUSLITE_230_CLOSEOUT_2026-10-19.md).
The old rc1 integration worktree is a retained historical candidate; its
scientific changes are absorbed and it is not a pending main-branch merge.

Release 2.2.9 includes the following changes, absent from 2.2.8.
Historical local acceptance (2026-10-15): the published capability
rows below do not certify every boundary input. The current main branch
contains fixed-cell/MIC/IDPP repairs, Relax input fallback and distribution notice
fixes, synchronized API progress failures, Dimer termination facts and explicit
stored-data TS validation. The CLI/API shared execution layer and backend provenance are integrated,
with 1284 default tests, 27 real-MPI checks, isolated distribution inspection
and clean-wheel public API checks passing. Validation distinguishes known
vibration scope from unknown evidence; no calculation automatically certifies a TS.
See the [architecture and distribution review](ARCHITECTURE_DISTRIBUTION_REVIEW_2026-10-13.md)
for acceptance evidence and validation limits. No new package release,
backend/GPU run or platform acceptance is claimed.

The nine issues found in the independent review of `5ad5de5` (one P1 and eight
P2) are repaired and fast-forwarded into local `main` as `9e5f3b8`
(included in the 2.2.9), from `fix/review-nine-findings-20261016`: segmented NEB restart, IRC default output
and frozen results, DMF endpoint element order, manifest path bases, bench
directory/resource/case-ID boundaries, and constant-potential facts-v1 compatibility.
The final source passes 1347 default tests (32 skipped), 27 opt-in local MPI
checks, and clean-wheel public API gates; independent source review passed.
See [the repair acceptance report](MAIN_REVIEW_FIXES_2026-10-16.md) for evidence
and migration limits. The original findings remain in
[the historical review](MAIN_BRANCH_REVIEW_2026-10-15.md#6-当前实现独立系统复核2026-10-15).
These fixes are not part of the published 2.2.8 package; “Supported” below
continues to describe capability availability.

The [2026-10-16 architecture review](ARCHITECTURE_REVIEW_2026-10-16.md) confirms
the current modular package, shared scientific execution and runtime/bench split
as the basis for continued development. Its recommended local contract repairs
are recorded in the acceptance report above; this does not certify a release or
close platform P6.

The 2.2.7 package was published from tag `v2.2.7` at commit
`af9c8fa6b82386a52055a5b3dced82eb3d8f4ef2`. It adds shared-GPU-node resource binding and run evidence (the
`runtime` section, worker isolation, `runtime_evidence.json`), the benchmark
toolchain (`atst_tools.bench`), the ABACUS `omp`/manifest thread-budget fixes and
the constant-potential development candidate. Local validation, the
documentation gate, the clean-wheel API gate and the GitHub release preflight
passed; platform runtime evidence remains pending. The 2.2.6 package was published from
tag `v2.2.6` at commit
`a633f06b9375e3f34dbd87f90f40ec6271b29cd4`. It records optimizer convergence
facts across workflows, writes owner manifests for Sella/Relax/AutoNEB, and
unifies runtime diagnostics to English. GitHub Tests, abacuslite, PyPI
publication and the official no-cache clean-install verification passed;
SIF/SAI/platform runtime evidence has not been performed. The 2.2.5 release
evidence stays in the 2.2.5 release notes.

Development history (2026-09-20; subsequently included in 2.2.7): CCQN records effective bond endpoints,
direction provenance and initial-structure identity. PRFO evaluates the previous
step with its original Hessian and applies the updated radius to the next step.
Local unit validation: 875 passed, 2 skipped; no DFT or deployment claim.
See the CCQN section of [CONFIG_REFERENCE](../user/CONFIG_REFERENCE.md).

Development history (2026-09-20; candidate code subsequently shipped in 2.2.7/2.2.8): the constant-potential candidate adds an
ABACUS calculator decorator, a fixed-geometry single-point/serial-scan workflow,
and fixed-cell `relax`/`neb` routing for the explicit `compensated_gate`
boundary. The legacy `reference_fcp` boundary remains reference-only. Local
ATST behavior checks, bounded numerical profile evidence and a 2026-09-21 SAI 4V100
joint run with the `runtime` section exist (gate single point plus three-point scan,
2/2 cases). It has been included since 2.2.8 and remains a development candidate
in 2.3.1, with no Paimon/public tool-chain or deployment acceptance claim. See the
[CP configuration section](../user/CONFIG_REFERENCE.md#211-constant-potential-evaluation-development-candidate)
for fields, artifacts, and failure semantics.

| Feature | Status | Description | Notes |
| :--- | :--- | :--- | :--- |
| **Relax** | ✅ Supported | Geometry Optimization | Uses ASE optimizers. |
| **Constant-potential** | 🧪 Development candidate | ABACUS-only fixed-electron electronic-number loop and fixed-geometry serial scan; fixed-cell `relax`/`neb` may use the explicit compensated-gate boundary. | `reference_fcp` is reference-only; compensated targets use custom `mu` in eV and same-run gate/dipole density facts. No cell-relax/stress, AutoNEB, D2S, MD, two-Fermi, or explicit `nupdown` scope. It remains a development candidate in 2.3.1; the public Paimon chain remains pending. The 2.2.9 CP-01 repair restores the legacy `profile_status` label on facts-v1 output and tests actual 2.2.8 readers; transition facts without the label require validated republishing for old consumers. The 2026-10-20 Pt/LCAO single-point `reference_fcp` parity on ABACUS LTS 3.10.1 passed with abacuslite 1.2.1 (see [the parity report](ATST_CP_PARITY_2026-10-20.md)); 1.2.1 includes the entropy-aware `free_energy` fix. The same fixture also completed both `reference_fcp` and `compensated_gate` on develop/latest-format builds `v3.11.0-beta6` (job 1864431), `v3.11.0-beta9` (job 1864468), and a locally built `origin/develop` `7be5df975` / `v3.11.0-beta11` (job 1864773); these are execution/format checks, not FCP v2 parity or scientific acceptance. |
| **Vibration** | ✅ Supported | Frequency Analysis and TS validation | Finite difference method with JSON results, TS validation, and artifact manifest support. |
| **NEB** | ✅ Supported | Nudged Elastic Band | CI-NEB, two-stage NEB, endpoint single-point repair, optional endpoint relaxation, native ASE selector, artifact manifest, ABACUS STRU inputs for `atst neb make`, and MPI image-level parallelism are supported. |
| **AutoNEB** | ✅ Supported | Automated NEB | Adaptive image handling, native ASE selector, endpoint single-point repair, and MPI image-level parallelism are supported. |
| **Dimer** | ✅ Supported | TS Search | Min-mode following. |
| **Sella** | ✅ Supported | Saddle Point Finder | Robust optimization. |
| **CCQN** | ✅ Supported | Constrained Cone Quasi-Newton TS Search | Standalone workflow and D2S refinement option, including reactive-mode enumeration, product alignment, diagnostics, mode manifest, and artifact manifest. |
| **D2S** | ✅ Supported | Double-Ended to Single | `atst run` dispatches endpoint optimization, rough DyNEB, then Dimer, Sella, or CCQN refinement. Experimental `rough_method: dmf` can replace rough DyNEB but is not a supported production default. |
| **IRC** | ✅ Supported | Intrinsic Reaction Coordinate | Sella backend and descent backend are supported, with controlled boundary diagnostics and artifact manifests. |
| **MD** | ✅ Supported | Molecular Dynamics | Supports ASE-driven MD with ABACUS/DP calculators and ABACUS-native MD input/run/output orchestration. Native MD selects the exact `OUT.<suffix>` from its INPUT, retains flat-output compatibility, and fails rather than reporting completion when no frames are produced. Bounded current-develop native-MD handoff passed on SAI job 1820042, including physical cell, periodic boundaries and raw trajectory checks; this does not establish scientific convergence or indexed-package release. |
| **Independent abacuslite required backend** | ✅ Supported (2.3.1) | One installed backend for calculator, ABACUS input/output, reverse configuration, native MD and gate compensation | Requires published `abacuslite>=1.2.1,<2`; missing or broken imports fail closed. Historical `atst_tools.external.ASE_interface.abacuslite` paths forward to the same installed modules for Paimon and other callers. The installed-wheel external consumer lane runs in CI, with a weekly schedule and an optional abacuslite ref for coordinated validation. Prior paired vendored/external wheel evidence remains frozen audit history. No new scientific convergence or TS claim follows from this migration. |
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
| **Runtime Device Binding & Evidence** | ✅ Supported | Optional device selection, thread budgets and run evidence | The optional `runtime` YAML section plus `--devices/--binding/--threads/--telemetry` on `atst run` and the API runner: requests are resolved against the caller's binding and trusted allocation facts (`ATST_ALLOCATION_DEVICES`), `round_robin` is fail-closed for unverified pools, and `runtime_evidence.json` records the environment triple, device facts, process-scope counters, host samples and the sampled memory peak. Legacy calls without runtime keys are unchanged. Local real-GPU evidence exists, and the SAI V100 P5 matrix completed on 2026-09-21 (DP matrix 12/12, ABACUS examples 4/4, multi-card NEB, 8 ranks on one card, 8 graphs x 8 GPU cards, and the host/SIF pairing; see sections 5h/5j of the SAI validation report). The DP TF artifact dimension is not specifically supported and its dedicated validation is not planned. The site `mpiexec` + runtime rebind fix passed three consecutive site runs (see the archived MPI repeat evidence); remaining site/upstream coordination is distinct from the reproduced runtime fix. P6 platform migration and downstream acceptance remain open. The 2026-10-17 parent candidate reached dev, but real entry checks stopped at missing SIF dependency APIs and unsupported Lustre no-replace publication; no GPU invocation was submitted. Flow registry returned HTTP 500, and prod remains unpublished (P6 takeover report, section 8). On 2026-10-16, the parent repository isolated branch fixed serial-image inner-MPI allocation/mapping handoff and TS-guess image counting (182 focused checks and independent review); this is local implementation evidence, not GPU/platform acceptance (see the P6 takeover report, section 7). Read-only model sharing and restart regression were closed on 2026-09-22 with a read-only model directory, a shared-worker DP batch, per-attempt attribution and a same-workdir restart on SAI V100; on the same site MPS is enabled, so per-process attribution stays `unavailable` by design instead of being inferred (see the P2 closeout report). A 2026-09-21 joint run with the constant-potential workflow passed on SAI; it also exposed a schema default that kept `runtime.threads` from reaching the ABACUS backend, which was fixed in `5e26789` before 2.2.7 and re-verified on site (SAI job 1436782; tracked in the GPU tuning review map). New in 2.2.7; additive. |
| **atst prepare** | ✅ Supported | Reverse config generation from an ABACUS run directory | `atst prepare` and `atst_tools.api.build_config_from_abacus_dir` generate runnable transition YAML (currently NEB) from `INPUT`/`STRU`/`KPT`/`PP`/`ORB`. User-controlled ABACUS values are kept verbatim with three technical floors (`calculation`→`scf`, `cal_force`→`1`, KPT line-mode rejected) and an optional endpoint energy+forces gate (`--no-gate` skips it). New in 2.2.3; backward-compatible. |
| **Trajectory Stress Retention** | ✅ Supported | NEB/AutoNEB trajectories carry per-image stress when the ABACUS calculator requests it | `cal_stress=1` retains per-image stress on NEB/AutoNEB trajectories (serial and image-parallel); `cal_stress=0` behavior is identical to 2.2.x. New in 2.2.3; additive. |
| **GA** | ❌ Not Supported | Genetic Algorithm | ASE 3.28.0 moved GA implementation to the standalone `ase-ga` project; ATST-Tools does not expose GA workflows. |
