# ATST-Tools 2.3.2 Release Notes

**Version**: 2.3.2
**Date**: 2026-10-20
**Status**: Published on PyPI and GitHub
**Source**: `main`, release commit `05773588e8722a7d45c249cf32480f7ece8e2066`
**Tag**: `v2.3.2`

## Summary

This patch release adds experimental fixed-cell constant-potential molecular
dynamics to the existing ASE MD workflow and compensated-gate calculator.

- Use `calculation.type: md`, `driver: ase` and the existing
  `calculator.constant_potential` mapping with one scalar custom target.
- Reuse ASE VelocityVerlet, Bussi, Langevin and NVTBerendsen. Converge the
  electronic state and electron number before returning each nuclear force.
- Save complete-step checkpoints, electron/potential diagnostics and grand
  energies. Validate restart identity and output tails before continuing.
- Read saved CP trajectory facts with offline MD postprocessing.
- Add a Pt example materializer, analytical regressions and real-GPU validation.
- Include the English documentation structure and Chinese onboarding guide
  maintained on main after 2.3.1.

See the [configuration reference](../user/CONFIG_REFERENCE.md#experimental-constant-potential-md)
and [Pt example](../../examples/19_constant_potential_Pt/README.md).

## Compatibility

- Package version: `2.3.2`.

This is an additive release. Existing YAML, CLI and public Python API entry
points remain available; the CP MD path composes existing MD and calculator
settings. The backend requirement remains `abacuslite>=1.2.1,<2`. No abacuslite
source change or unreleased ABACUS socket extension is required.

CP MD accepts only fixed-cell ASE dynamics with `compensated_gate`, a custom
reference and one target. Native ABACUS CP MD, NPT/variable cell, target lists
and `reference_fcp` MD remain unsupported. The `nve` selector means
unthermostatted nuclear VelocityVerlet; electron exchange makes it different
from a fixed-electron-number microcanonical calculation.

Restart continues a saved physical state. Thermostat RNG state is not saved,
so stochastic continuation is not guaranteed to reproduce an uninterrupted
trajectory. Inconsistent output tails are rejected rather than repaired.

## Validation

Implementation commit `156d0c5` passed independent review after repairs to
native-driver rejection and restart nuclear-state/energy validation.

- Full offline source suite: **1450 passed, 37 skipped**, zero failures.
- Clean installed-wheel external consumer: **1452 passed, 35 skipped**, zero
  failures or deselections; installed-origin and report verification passed.
- Wheel/sdist build, strict Twine and clean-wheel public API checks passed.
  These pre-release checks used the implementation before the version bump;
  release CI must verify the final 2.3.2 commit and artifacts.
- Analytical tests exercise the electronic solve, minimized grand-potential
  forces, timestep refinement, integrators, constraints and deterministic
  restart. Real ABACUS LTS 3.10.1 / abacuslite 1.2.1 / V100 runs verify short
  NVT/Verlet trajectories and split continuation.
- Final-source GPU job `1873489` completed successfully. The maximum CP
  residual was below 2.46e-6 eV against a 1e-5 eV tolerance; split versus
  continuous positions differed by at most 7.41e-12 Angstrom.

The [dated validation report](../reports/ATST_CP_MD_VALIDATION_2026-10-20.md)
records exact source, environments, jobs and scientific limits. The real
backend runs are short: they do not establish long-time sampling or clean
second-order timestep convergence. Constant-potential workflows retain their
experimental status.

## Publication

- [PR #47](https://github.com/QuantumMisaka/atst-tools/pull/47) merged after all
  six required checks passed in [candidate CI](https://github.com/QuantumMisaka/atst-tools/actions/runs/39392303090).
  The merged tree equals reviewed candidate `ab45fc6`.
- [Main CI](https://github.com/QuantumMisaka/atst-tools/actions/runs/39392582246)
  passed on release commit `0577358`. The local exact-tag readiness check passed.
- [Tagged publication](https://github.com/QuantumMisaka/atst-tools/actions/runs/39392692089)
  passed readiness, full tests, documentation governance, build, strict Twine
  and clean-wheel API checks. PyPI upload completed after the existing
  five-minute environment timer.
- [PyPI 2.3.2](https://pypi.org/project/atst-tools/2.3.2/) was uploaded on
  2026-10-20 at 08:27 UTC (16:27 UTC+08:00). Indexed artifact sizes and SHA256
  digests match the downloaded tag artifacts. The
  [GitHub release](https://github.com/QuantumMisaka/atst-tools/releases/tag/v2.3.2)
  carries those same files.
- The tag wheel's 93 Python files match the release source. With abacuslite
  1.2.1 and ASE 3.28.0 in an isolated Python 3.10 environment, its full
  installed-consumer suite collected 1487 tests: **1452 passed, 35 skipped**,
  zero failures or deselections. Installed-origin/report verification passed.
- A separate no-cache PyPI install into a new Python 3.10 environment passed
  `pip check`, CLI/API version checks (`2.3.2`), external-backend identity and
  comparison of all 93 installed Python files against the tag wheel. This
  unpinned dependency resolution selected abacuslite 1.2.2 and ASE 3.29.0;
  the lower-bound consumer and GPU evidence above used abacuslite 1.2.1.

| Artifact | Bytes | SHA256 |
| --- | ---: | --- |
| `atst_tools-2.3.2-py3-none-any.whl` | 362240 | `25eaa0f4619d96c1d8540a9cc266196a1b0816feb3d3063c34cda1a5a7f32a83` |
| `atst_tools-2.3.2.tar.gz` | 321117 | `d5565bf75becba3757422ca58e8c4897e0c48e1fd30ac38a12a2a631160d7f15` |

Raw logs, tag artifacts, consumer report and PyPI metadata are retained as
`release232-*` under `/home/james/scratch/atst-abacuslite-release-20261020/`.
The original checkout and scientific evidence remain in place. This release
does not update downstream Toolbox deployments. Publication records added to
main after the tag do not alter the tagged artifacts.
