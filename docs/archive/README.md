# Archive

This directory keeps historical material that is no longer current guidance.
Active documentation lives in `docs/index.md`, `docs/user/`,
`docs/developer/`, `docs/releases/` and `docs/reports/`.

## Deletion policy

A full document may be deleted when all of the following hold:

- no active Markdown link points to it;
- it has no unique validation, environment, scientific or release evidence;
- its conclusions are absorbed into maintained docs, tests or release notes;
- its identity and disposition are recorded in this index.

Git history remains the full-text archive after deletion. A one-line summary
here is the audit pointer.

## Consolidated governance records

The former `DOCS_ARCHITECTURE.md` and `DOCUMENTATION_STANDARDS.md` were
consolidated into [the documentation guide](../developer/DOCUMENTATION.md).
Their full historical copies were removed after this index was written.

The former documentation status ledger was condensed into
[HISTORY.md](HISTORY.md). Release notes and validation reports remain the
authoritative evidence for their dates.

## Release notes

Release notes for `2.0.0` through `2.2.6` are in [releases/](releases/).
`2.2.7` through `2.2.9` and the current `2.3.0` candidate remain in
`docs/releases/`.

## Reports

Reports dated before 2026-09 were moved to [reports/](reports/). Current
evidence remains in `docs/reports/`; the feature matrix is the status entry.

## Pending-delete history

The following files were reviewed as superseded and removed from active
documentation. Full text remains in Git history. The retained conclusions are
covered by the current implementation, tests, release notes or maintained
references named in the former ledger.

- `plans/2026-08-04-abacuslite-force-read-and-cross-repo-plan.md` — delivered;
  cross-repository follow-up remains with the owning repositories.
- `plans/2026-09-17-workflow-convergence-handoff-plan.md` — delivered in 2.2.6;
  handoff contracts live in the API reference and convergence fixtures.
- `specs/2026-08-12-atst-prepare-reverse-and-trajectory-stress-design.md` —
  delivered in 2.2.3; current behavior is covered by the CLI/API reference.
- `specs/2026-09-04-user-distribution-and-mpi-installation-design.html` —
  delivered in 2.2.4; installation boundaries live in the user guide.
- `specs/2026-08-04-abacuslite-cross-repo-and-force-read-design-review.md` —
  conclusions absorbed into the maintained abacuslite design.
- `plans/2026-10-16-review-findings-remediation.md` — nine review findings
  repaired locally; evidence retained in `MAIN_REVIEW_FIXES_2026-10-16.md`.
- `plans/2026-10-13-ase-ts-dual-interface-remediation.md` — local remediation
  completed; evidence retained in the system review.
- `plans/2026-07-02-abacuslite-backend-upstream-issue-fixes.md` — implemented;
  upstream synchronization is tracked by `PATCHES.md` and the snapshot baseline.
- `plans/2026-07-05-atst-cli-banner.md` — shipped with CLI tests and user docs.
- `plans/2026-07-05-ci-test-development.md` — CI contracts live in the
  workflow and handover checklist.
- `plans/2026-07-22-abacus-agent-atst-api-adoption.md` — runner adoption shipped.
- `plans/2026-09-04-user-distribution-and-mpi-installation.md` — delivered in
  2.2.4; conclusions live in README, user guide and release notes.
- `plans/2026-09-04-github-governance-and-release-gates.md` — repository-side
  gates implemented; administrator settings remain external.
- `plans/2026-09-16-abacuslite-periodic-frame-validation-225-release.md` —
  delivered in 2.2.5; evidence in release notes.
- `plans/DMF-integrate-and-test-plan.md` — DMF shipped as experimental; current
  boundaries are in the feature matrix.
- `plans/Direct-MaxFlux-方法调研.md` — temporary note superseded by the DMF
  research report.
- `plans/native-ase-backend.md` — implemented; current behavior is covered by
  the backend selector and tests.
- `reports/PROJECT_REFACTOR_REVIEW_2026-05-15.md` — old refactor baseline.
- `reports/USER_EXPERIENCE_REINFORCEMENT_2026-05-15.md` — UX work absorbed.
- `reports/CY_PT_AUTONEB_MAIN_REPRODUCTION_REVIEW_2026-05-18.md` — superseded
  by the Issue #25 fix validation.
- `reports/CY_PT_AUTONEB_MAIN_ALIGNED_LTS3101_VALIDATION_2026-05-19.md` —
  early unsuccessful validation, superseded.
- `reports/CY_PT_AUTONEB_FAILURE_ROOT_CAUSE_REVIEW_2026-05-18.html` — root
  cause absorbed by the final fix report.
- `reports/ISSUE_25_AUTONEB_FMAX_REVIEW_2026-05-18.md` — pre-fix assessment,
  superseded by the final fix report.
- `reports/ISSUE_25_AUTONEB_SAI_VALIDATION_2026-05-18.md` — early validation,
  superseded.
- `reports/ZN_MIGRATION_NEB_ABACUS_VALIDATION_2026-05-26.md` — initial Zn
  validation, superseded by the segmented runtime report.
- `reports/EXAMPLES_INITIAL_GUESS_AUDIT_2026-05-26.md` — absorbed by perturbed
  input validation and reference results.
- `reports/UNIT_TEST_COVERAGE_REVIEW_2026-05-25.md` — coverage review now lags
  later P0/P1, MPI and DP test additions.
- `trae/documents/CCQN-plan.md` — temporary CCQN plan, absorbed.
- `trae/documents/NEB_para_abacuslite.md` — temporary merge note, superseded.
- `trae/documents/NEB_parallel_imple.md` — temporary MPI implementation plan.
- `trae/zn-neb/README.md`, `Zn-NEB.md`, `checklist.md`, `spec.md`, `tasks.md` —
  superseded Zn planning material.

## Completed CP MD implementation

The [2026-10-20 implementation plan](developer/2026-10-20-constant-potential-md.md)
records the unreleased fixed-cell ASE CP MD delivery. Lasting contracts are in
the configuration reference and the dated MD validation report.
