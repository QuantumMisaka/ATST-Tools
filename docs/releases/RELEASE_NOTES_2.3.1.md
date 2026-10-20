# ATST-Tools 2.3.1 Release Notes

**Version**: 2.3.1
**Date**: 2026-10-20
**Status**: Release candidate prepared for PyPI publication
**Branch**: `feat/cp-reference-parity` (pending merge to `main`)
**Tag**: `v2.3.1` (assigned at release time)

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

- PyPI publication and artifact hashes are recorded after the tagged
  Package/Release workflow completes.
