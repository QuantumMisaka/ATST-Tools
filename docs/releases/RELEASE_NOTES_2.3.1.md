# ATST-Tools 2.3.1 Release Notes

**Version**: 2.3.1
**Date**: 2026-10-20
**Status**: Release candidate prepared for PyPI publication
**Branch**: `main`
**Tag**: `v2.3.1` (assigned at release time)

## Summary

This patch release corrects the constant-potential energy handoff against the
external FCP v2 reference implementation and adds reproducible parity tooling.

- Require `abacuslite>=1.2.1,<2`; 1.2.1 exposes the entropy-corrected
  `free_energy` needed by the CP grand potential.
- Use the force-consistent energy for the CP `energy` result and `E_KohnSham`
  for the CP `free_energy` result, matching FCP v2's result-key convention.
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
the entropy-aware `free_energy` handoff. No YAML schema field is removed or
renamed. `BACKEND_SOURCE` remains always `external`.

On smeared systems, CP `energy` and `free_energy` values change to the corrected
FCP v2 convention; this is a bug fix, not a new scientific calibration.
`constant_potential` remains a development candidate.

## Validation

- Default unit suite completed without failures.
- Documentation governance passed.
- SAI LTS job `1864305` passed FCP v2 single-point `reference_fcp` parity with
  `ok: true` and no hard mismatches.
- SAI develop jobs `1864431` (`v3.11.0-beta6`), `1864468` (`v3.11.0-beta9`) and
  `1864773` (locally built `origin/develop` `7be5df975`, `v3.11.0-beta11`)
  completed both `reference_fcp` and `compensated_gate`.

## Publication

- PyPI publication and artifact hashes are recorded after the tagged
  Package/Release workflow completes.
