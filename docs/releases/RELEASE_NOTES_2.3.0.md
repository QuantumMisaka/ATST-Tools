# ATST-Tools 2.3.0 Release Notes

**Version**: 2.3.0
**Date**: 2026-10-19
**Status**: Source candidate (unpublished)
**Branch**: `main`
**Tag**: Pending

## Summary

This minor release makes the independently maintained
[`abacuslite`](https://github.com/QuantumMisaka/abacuslite) package the sole
ABACUS ASE backend. ATST-Tools keeps its workflow, YAML and CLI responsibilities;
calculator construction, ABACUS input/output and parser behavior move to the
installed package.

- Declare `abacuslite>=1.1.0,<2` as a required runtime dependency.
- Remove vendored selection: missing or broken installations now fail with a
  clear diagnostic instead of falling back to a bundled copy.
- Route calculator, INPUT/STRU/KPT I/O, output readers, reverse configuration,
  native MD and gate compensation through the installed package.
- Keep `atst_tools.external.ASE_interface.abacuslite` as thin forwarding paths
  for Paimon and other existing callers.
- Add an installed-wheel external consumer lane, a weekly dependency check and
  an optional abacuslite ref for coordinated cross-repository validation.
- Retain the configured ABACUS executable, including `env` and MPI launcher
  wrappers, for version probes.

## Compatibility

- Package version: `2.3.0`.

Existing CLI, YAML and public Python API entry points remain available. The
mandatory `abacuslite` dependency is the deliberate compatibility boundary:
offline or air-gapped installations that relied on the bundled snapshot must
now provide `abacuslite>=1.1.0,<2`. `BACKEND_SOURCE` is always `external`;
callers that branched on the former `vendored` value must adapt. Historical
import paths continue to resolve to the installed modules. No YAML schema field
is removed or renamed.

## Validation

The source candidate passes the default unit suite with five environment-gated
skips, documentation governance, and the installed-wheel external consumer
suite. Remote GitHub tests and the publication gate are recorded separately
after merge. Publication, tag creation and PyPI upload are not implied by this
source release note.

## Related documents

- [ABACUSLite wrapper guide](../user/ABACUSLITE_WRAPPER_GUIDE.md)
- [Feature status matrix](../reports/FEATURE_STATUS_MATRIX.md)
- [2.2.9 release notes](RELEASE_NOTES_2.2.9.md)
