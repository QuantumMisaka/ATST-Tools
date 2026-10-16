# ATST-Tools Documentation

This is the maintained entry point for ATST-Tools documentation. Files under
`.trae/` are temporary working notes and are not part of the formal
documentation set.

## Start Here

| Your starting point | Recommended route |
| --- | --- |
| New user without a calculator runtime | [Quick start](user/USER_GUIDE_CN.md#3-10-分钟快速上手): obtain matching examples, validate YAML, then create/read local files. |
| Ready to run ABACUS or DP | [Examples](../examples/README.md#usage): check resources and data, enter the case directory, run, then inspect convergence. |
| Agent or Python integration | [Python API](user/PYTHON_API_REFERENCE.md): use typed results or the process runner and keep each calculation in its own directory. |
| New contributor | [Developer setup](developer/HANDOVER.md#11-新-checkout-与-cpu-开发准备), then [source architecture](../src/atst_tools/README.md). |
| Project maintainer exploring implementation | [Source architecture](../src/atst_tools/README.md) for package responsibilities, execution relationships, and the purpose of runtime and bench. |

The current release is 2.2.8. Sections marked **unreleased development tree**
describe source changes that are not available by installing that release.
Package versions do not establish scientific convergence or platform acceptance.

## User Path

Use this path when you want to install ATST-Tools, choose an example, validate a
YAML file, and run a workflow.

从 PyPI 阅读 README 时，请使用 README 顶部的 GitHub/Gitee 绝对用户文档入口；仓库内浏览继续使用下列本地链接。

1. [README quick start](../README.md) for project scope, supported workflows,
   installation, and minimal commands.
2. [Chinese user guide](user/USER_GUIDE_CN.md) for a 10-minute onboarding path
   covering local file operations before a separately configured backend run.
3. [Examples overview](../examples/README.md) for runnable NEB, AutoNEB, Dimer,
   Sella, CCQN, D2S, Relax, Vibration, IRC, MD, experimental DMF, DP, and MPI
   examples, plus the constant-potential development validation fixtures.
4. [CLI reference](user/CLI_REFERENCE.md) for `atst run`, `atst config`,
   `atst abacus`, and lightweight helper commands.
5. [Stable Python API reference](user/PYTHON_API_REFERENCE.md) when embedding
   the same schema-backed workflows or CCQN in another Python program, or when
   an external host needs the installed API runner and JSON handoff protocol.
6. [Configuration reference](user/CONFIG_REFERENCE.md) for hand-written YAML
   semantics and common patterns.
7. [YAML input variables](user/YAML_INPUT_VARIABLES.md) for the generated
   schema field table.
8. [ABACUSLite wrapper guide](user/ABACUSLITE_WRAPPER_GUIDE.md) for ABACUS
   backend boundaries and image-level MPI configuration.
9. [Feature status matrix](reports/FEATURE_STATUS_MATRIX.md) for supported,
   partial, experimental, and unreleased workflow boundaries.

## Developer Path

Use this path when you need to modify a workflow, YAML schema, calculator
backend, CLI command, example, or documentation.

1. [Source architecture guide](../src/atst_tools/README.md) for package
   responsibilities, execution relationships, and change entry points.
2. [Handover checklist](developer/HANDOVER.md) for routine maintenance tasks.
3. [Example validation operations](developer/EXAMPLE_VALIDATION_OPERATIONS.md)
   for maintainer checks, site execution, curated-output provenance, and
   evidence reports.
4. [YAML input governance](developer/YAML_INPUT_GOVERNANCE.md) for schema,
   generated docs, examples, and config tests.
5. [Documentation standards](developer/DOCUMENTATION_STANDARDS.md) for metadata,
   lifecycle, report levels, archive rules, and verification commands.
6. [Documentation architecture](developer/DOCS_ARCHITECTURE.md) for directory
   responsibilities and target audiences.
7. [Maintained atst-cli skill](skills/atst-cli/SKILL.md) for operational CLI
   usage and validation snippets.
8. [PyPI release automation](developer/PYPI_RELEASE_AUTOMATION.md) for release
   publishing and the local readiness gate.
9. [Governance and release gates](developer/GOVERNANCE_AND_RELEASE_GATES.md)
   for mechanical checks, cross-family review, and administrator/post-push
   responsibilities.

## Project Manager Path

Use this path when you need to judge current feature support, validation
evidence, documentation health, release scope, or cleanup priorities.

For a map of what the implementation contains and why each directory exists,
start with the [source architecture guide](../src/atst_tools/README.md).

1. [Feature status matrix](reports/FEATURE_STATUS_MATRIX.md) for current
   supported, partial, and unsupported capabilities.
2. [Documentation governance report](reports/DOCUMENTATION_STATUS_REPORT.md) for
   the active documentation ledger, report levels, archive state, and
   pending-delete status.
3. [2.2.8 release notes](releases/RELEASE_NOTES_2.2.8.md) for the current
   version-level delivery scope and recorded tag, CI, and PyPI evidence;
   official no-cache clean-install, CLI/API, and dependency verification passed.
   The previously published [2.2.5 release notes](releases/RELEASE_NOTES_2.2.5.md)
   and [2.2.4 release notes](releases/RELEASE_NOTES_2.2.4.md) remain available
   for historical publication evidence.
4. Current validation reports linked from the governance report when judging a
   specific feature, backend, or environment.

## Backend Policy

ATST-Tools consumes the official ABACUS ASE interface through `abacuslite`. The
resolver first tries an independently installed `abacuslite` package, then
falls back to the vendored snapshot under
`src/atst_tools/external/ASE_interface/abacuslite`. The repository ships the vendored snapshot so a separate installation is not
required. An independently installed package takes precedence and may behave
differently; the selected source should be considered when reproducing a run.

Maintenance is *this-repository-first*: fixes and hardening land in the
vendored snapshot first, are recorded in
`src/atst_tools/external/ASE_interface/PATCHES.md`, and are periodically synced
back upstream to `abacus-develop`. Upstream updates are pulled on a regular
cadence against the baseline tracked in
`src/atst_tools/external/ASE_interface/ABACUSLITE_SNAPSHOT.md`.

The long-term goal remains unchanged: once `abacuslite` is independently
installable with a stable release channel, ATST-Tools should move it to an
optional dependency or extra and retire the vendored fallback.
