# ATST-Tools 2.2.9 Release Notes

**Version**: 2.2.9
**Date**: 2026-10-17
**Status**: Published (official PyPI clean-install verified)
**Branch**: `main`
**Tag**: `v2.2.9` -> `ec75edf`

## Summary

This patch consolidates the scientific-contract and delivery repairs completed
after 2.2.8. CLI and Python API share scientific execution; backend provenance,
explicit stored-result TS validation and Dimer termination facts are retained.

- Repair fixed-cell/PBC preparation, periodic MIC/IDPP behavior, Relax input
  fallback, distribution notices and synchronized MPI progress failures.
- Restore complete segmented NEB restart bands; retain IRC default output names
  and independent frozen energy/force results; reject mismatched DMF endpoints.
- Record explicit manifest path bases for nested/external layouts.
- Reject conflicting benchmark work directories, resource overrides and unsafe
  case IDs; preserve configuration identity across worker working directories.
- Restore the legacy facts-v1 label required by 2.2.8 constant-potential readers.

## Compatibility

- Package version: `2.2.9`.

Existing CLI/YAML/public API entry points remain available. New facts and
manifest metadata are additive; malformed or ambiguous inputs now fail earlier.
Fixed-cell matching uses absolute `1e-5 Å` tolerances (`rtol=0`), so former
`1e-5`–`1e-4 Å` matches can be rejected. External abacuslite is not patched.
Old nested manifests with workdir-relative paths require an explicit path base
from their owner; default root layouts remain compatible. Consumers of nested
layouts must understand the new `path_base` field.

A trajectory already continued from a wrongly sliced segmented NEB band is not
automatically repaired. Restore its original complete band and matching make
settings. Facts-v1 transition frames without the old label must be validated and
republished before handing them to an unchanged 2.2.8 reader.

Constant-potential remains a development candidate and DMF experimental.
No workflow automatically certifies a transition state. Platform P6 acceptance
is tracked separately from package publication.

## Validation

The 2026-10-17 candidate passed the full default suite (1347 passed,
32 skipped), all 27 opt-in local MPI checks, distribution inspection, strict
Twine checks, documentation governance and clean-wheel public API verification.
Independent review confirmed the repair evidence matches the current source
(the package-version change is the only changed build input).
GitHub publish run `38994894225` completed successfully. The Python 3.12 matrix exposed an overly broad
IDPP test warning assertion under NumPy 2.5. A post-tag test-only correction
excludes upstream deprecation warnings while retaining scientific warnings;
the full IDPP file passes in Python 3.12. Production source is unchanged.
Official PyPI no-cache wheel download and clean-environment API/CLI/runner
verification passed. Post-tag Tests `38995247103` and abacuslite
`38995247138` passed. The tag remains immutable at `ec75edf`.

## Published artifacts

- Wheel: `atst_tools-2.2.9-py3-none-any.whl`, 415854 bytes, SHA-256
  `a62c8d346206226650ad4ced1a165e9a0a141b5ad4a5066e119c35d463d96a15`.
- Sdist: `atst_tools-2.2.9.tar.gz`, 366100 bytes, SHA-256
  `63c6829b733e0563c26d01e0674afef91922a28d2027110fe02b0679c580e20c`.

See [repair acceptance](../reports/MAIN_REVIEW_FIXES_2026-10-16.md) and
[scientific-contract review](../reports/ARCHITECTURE_DISTRIBUTION_REVIEW_2026-10-13.md).
