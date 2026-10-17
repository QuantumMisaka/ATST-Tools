# ATST-Tools 2.2.9 Release Notes

**Version**: 2.2.9
**Date**: 2026-10-17
**Status**: Release candidate (publication pending)
**Branch**: `main`
**Tag**: `v2.2.9` (pending)

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
Publication and official PyPI installation checks remain pending.

See [repair acceptance](../reports/MAIN_REVIEW_FIXES_2026-10-16.md) and
[scientific-contract review](../reports/ARCHITECTURE_DISTRIBUTION_REVIEW_2026-10-13.md).
