# ATST-Tools 2.3.0 Release Notes

**Version**: 2.3.0
**Date**: 2026-10-19
**Status**: Published on PyPI (2026-10-19)
**Branch**: `main`
**Tag**: `v2.3.0` at `b8c3218` (tag object `04334a7`)

## Summary

This minor release makes the independently maintained
[`abacuslite`](https://github.com/QuantumMisaka/abacuslite) package the sole
ABACUS ASE backend. ATST-Tools keeps its workflow, YAML and CLI responsibilities;
calculator construction, ABACUS input/output and parser behavior move to the
installed package.

- Declare `abacuslite>=1.2.0,<2` as a required runtime dependency.
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
now provide `abacuslite>=1.2.0,<2`. `BACKEND_SOURCE` is always `external`;
callers that branched on the former `vendored` value must adapt. Historical
import paths continue to resolve to the installed modules. No YAML schema field
is removed or renamed.

The validated 2.3.0 consumer baseline is `abacuslite==1.2.0` (PyPI, 2026-10-19);
the earlier `1.1.0` candidate is outside the supported range and is retained
only as historical acceptance evidence.

## Validation

The 2.3.0 candidate passed the default unit suite with five environment-gated
skips, documentation governance, and the installed-wheel external consumer
suite. Hosted tests and the publication gate are recorded under Publication
below.

The [integration closeout record](../reports/ABACUSLITE_230_CLOSEOUT_2026-10-19.md)
tracks the current distribution, installed-consumer and real-engine acceptance,
the published plugin baseline, and remaining delivery work.

The historical `abacuslite==1.1.0` candidate at `49c8626` passed installed-distribution acceptance
(1386 passed, 33 environment-gated skips) and SAI job 1848462 (18 stages,
13 real-engine cases, including parallel NEB and LTS 3.10.1). These bounded
execution checks do not establish scientific convergence. Exact wheel hashes
and retained evidence locations are recorded in the closeout report above.

The `abacuslite==1.2.0` baseline re-run reproduced installed-package results;
it did not rerun SAI job 1848462 or establish GPU acceptance for 1.2.0:
the installed-state suite collects 1419 tests with 0 failures and 35
environment-gated skips, and the installed-wheel external consumer suite
reports 1419 collected / 1386 passed / 33 skipped / 0 failed, identical to the
1.1.0 run. The 1.2.0 wheel SHA256 is
`3e5e245c0a03e8182c2c726d7ccad1c800c20250641abd4b81a8e22e892bc4a0`.

## Publication

- Tag `v2.3.0` points to `b8c3218ade7380116d3fbc2162cf11bb7a8dded8`
  (squash merge of #40); annotated tag object
  `04334a754b6d32a4b1caabc3a3cfd3c86d22dc12`.
- Publish workflow run
  [39196721131](https://github.com/QuantumMisaka/atst-tools/actions/runs/39196721131)
  (`push` on `v2.3.0`) completed successfully, so preflight and upload both
  passed.
- PyPI artifacts: `atst_tools-2.3.0-py3-none-any.whl` SHA256
  `80a5e2422176df48bfdcc1e88bd66a11cdcb953c735ae86fe62cb91004fd5729`
  and `atst_tools-2.3.0.tar.gz` SHA256
  `af317351ab95e2526301288619fc55fd8c316fbf6a3f9711b50b5fc75aef3bf2`
  (uploaded 2026-10-19T03:59:48Z and 03:59:50Z).
- Published metadata requires `abacuslite<2,>=1.2.0`.
- Fresh-venv check of the published package: version `2.3.0`,
  `BACKEND_SOURCE` reports `external`, and `pip check` is clean.

The published PyPI wheel was independently rechecked on 2026-10-19 with
`abacuslite==1.2.0`, Python 3.10 and ASE 3.28.0 in a fresh environment:
1419 collected / 1386 passed / 33 skipped / 0 failed / 0 deselected. Both
published distribution artifacts were supplied to the suite, and their SHA256
values match PyPI. The consumer report verifier, `pip check`, CLI version and
installed-package provenance checks passed. Raw evidence is retained under
`/home/james/scratch/atst-230-finalize-20261019-za20o8qf/`.

## Related documents

- [ABACUSLite wrapper guide](../user/ABACUSLITE_WRAPPER_GUIDE.md)
- [Feature status matrix](../reports/FEATURE_STATUS_MATRIX.md)
- [2.2.9 release notes](RELEASE_NOTES_2.2.9.md)
