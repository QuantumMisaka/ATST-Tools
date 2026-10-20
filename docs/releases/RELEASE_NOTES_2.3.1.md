# ATST-Tools 2.3.1 Release Notes

**Version**: 2.3.1
**Date**: 2026-10-20
**Status**: Published on PyPI and GitHub
**Source**: `main`, release commit `7e54df2c805b984f93b18546aa1c941b13e198b2`
**Tag**: `v2.3.1`

## Summary

This patch release corrects the constant-potential energy handoff against the
external FCP v2 reference implementation and adds reproducible parity tooling.

- Require `abacuslite>=1.2.1,<2` for the explicit `energy_without_entropy`
  handoff used by the reference-only FCP v2 comparison.
- In `reference_fcp`, use the no-entropy energy for CP `energy` and preserve
  force-consistent `free_energy`. Keep the backend energy relationship intact
  for `compensated_gate` and ordinary ASE calculations.
- Expose the evaluated CP capacitance diagnostic.
- Add the FCP v2 reference producer, normalized record extractor, comparator,
  source-reproduction test and opt-in real-engine parity test.
- Record LTS FCP v2 single-point parity and beta6/beta9/latest-develop CP smoke
  evidence.
- Streamline the README and align current user/developer documentation with the
  external-abacuslite architecture.

## Compatibility

- Package version: `2.3.1`.

Existing CLI, YAML and public Python API entry points remain available. The
mandatory `abacuslite` floor moves from `1.2.0` to `1.2.1` because 1.2.0 lacks
the explicit no-entropy energy handoff. No YAML schema field is removed or
renamed. `BACKEND_SOURCE` remains always `external`.

On smeared systems, `reference_fcp.energy` changes to the FCP v2 convention;
its `free_energy` and both `compensated_gate` energy keys retain their original
meaning. This is a bug fix, not a new scientific calibration.
`constant_potential` remains a development candidate.

## Validation

- Local source suite: 1399 passed, 37 environment/opt-in skips, zero failures
  (Python 3.10.19, ASE 3.28.0, corrected abacuslite 1.2.1 candidate wheel).
- Installed-wheel external consumer suite: 1436 collected, 1399 passed,
  37 skipped, zero failures/deselections; the consumer report verifier passed.
- Wheel/sdist build, strict Twine and clean-wheel public API gates passed.
  Evidence is retained in `/home/james/scratch/atst-abacuslite-release-20261020/`.
- Independent final review cleared both corrected code surfaces; the optional
  no-entropy fallback and its parity limit are documented in CONFIG_REFERENCE.
- Documentation governance passed.
- Historical SAI LTS job `1864305` passed FCP v2 single-point `reference_fcp` parity with
  `ok: true` and no hard mismatches.
- Historical SAI develop jobs `1864431` (`v3.11.0-beta6`), `1864468` (`v3.11.0-beta9`) and
  `1864773` (locally built `origin/develop` `7be5df975`, `v3.11.0-beta11`)
  completed both `reference_fcp` and `compensated_gate` on the pre-review
  candidate. These runs are not new GPU acceptance of the corrected release.
- Independent review caught the global free-energy mislabel and the
  compensated-gate key swap before integration. Regressions cover both
  boundaries; replay of retained LTS electronic facts checks the reference
  numerical convention without claiming a new engine run.

## Publication

- [PR #43](https://github.com/QuantumMisaka/atst-tools/pull/43) merged after
  all six jobs passed in [candidate CI](https://github.com/QuantumMisaka/atst-tools/actions/runs/39359925318).
  One initial dependency-install failure was caused by the newly published
  backend not yet appearing at that PyPI index endpoint; rerunning that job
  passed without a source or workflow change.
- The merged source tree equals reviewed candidate `6bcb195`. Local release
  readiness confirmed that `v2.3.1` points to the release commit above.
- [Tag publication run](https://github.com/QuantumMisaka/atst-tools/actions/runs/39360379146)
  passed readiness, full pytest, documentation governance, build, strict Twine
  and clean-wheel API checks, then uploaded the checked artifacts after the
  configured five-minute PyPI environment wait.
- [PyPI 2.3.1](https://pypi.org/project/atst-tools/2.3.1/) was uploaded on
  2026-10-20 at 02:50 UTC (10:50 UTC+08:00). Both indexed digests match the
  downloaded tag artifacts. The
  [GitHub release](https://github.com/QuantumMisaka/atst-tools/releases/tag/v2.3.1)
  carries those same two files.
- The tag wheel's 92 Python files match the reviewed source. Installed with
  the published abacuslite 1.2.1 in an isolated Python 3.10 / ASE 3.28.0
  environment, the full external consumer suite collected 1436 tests:
  **1402 passed, 34 environment-gated skips, zero failures/deselections**.
  This includes both distribution inspections and the pinned FCP v2 source
  reproduction test. The installed-origin/report verifier passed.
- A subsequent no-cache install of `atst-tools==2.3.1` from the PyPI index
  passed `pip check`, reports `atst 2.3.1`, and resolves both packages inside
  the isolated environment with `BACKEND_SOURCE=external`.
- Evidence: `published-consumer.json` and its pytest/attestation logs,
  `atst-pypi-install.log`, `atst-tag-artifacts.json`, `atst-pypi.json` and
  `atst-publish-logs.zip` in the evidence directory listed above.

| Artifact | Bytes | SHA256 |
| --- | ---: | --- |
| `atst_tools-2.3.1-py3-none-any.whl` | 357851 | `ae63a1f83d572eeec96f19f181d0f5521165645d4c2b4c677fbb6ad5206f899f` |
| `atst_tools-2.3.1.tar.gz` | 320773 | `282878603e9fc69af2aa5e523b74901cf35bd54344b7bb2c57eb8728497542d2` |

## Remaining independent work

ABACUS upstream PR #8175 remains open; its engine/harness refs must be rebound
and revalidated after upstream merge. GPU P6's same-input 1/4-card NEB and
remaining mixed-parallel/default-path acceptance stay open. The parent
Toolbox's 1.2.76 dev/prod GaN acceptance is recorded separately; this package
release does not update its dependency gitlinks, SIF or deployed runtime.
The original checkouts and retained scientific worktrees remain in place.

## Documentation follow-up (2026-10-20)

The repository README now focuses on capabilities, installation and a first
example. Detailed setup stays in the user references; current support and
historical/platform evidence are linked from the feature matrix. User/API
references now identify 2.3.1 and stored-data TS screening as released, while
retaining the experimental constant-potential and DMF boundaries. This is a
documentation-only follow-up on main; the published tag and artifacts are
unchanged, and no CLI, YAML or Python API behavior changes.

The primary language for active documentation is now English, including the
source architecture and developer checklists. The Chinese user guide remains
a complete onboarding path, with advanced details linked to the English
references. The language policy and historical-material exceptions are defined
in [Documentation](../developer/DOCUMENTATION.md#language-policy). English
setup and result-reading paths no longer depend on the Chinese guide.
