# Documentation

This is the single maintenance guide for the ATST-Tools documentation tree.

## Entry points

| Path | Purpose |
| --- | --- |
| `README.md` | Installation, quick start and product scope |
| `docs/index.md` | Start-here, user, developer and status navigation |
| `docs/user/` | User guides, CLI/YAML/API references |
| `docs/developer/` | Development, validation, release and maintenance guides |
| `docs/releases/` | Version-level release notes |
| `docs/reports/` | Current evidence and status; not a required reading path |
| `docs/superpowers/` | Design specs and active plans |
| `docs/archive/` | Historical material kept for audit |
| `docs/skills/` | Reusable operational references |

Validation and site operations are documented in
[EXAMPLE_VALIDATION_OPERATIONS.md](EXAMPLE_VALIDATION_OPERATIONS.md).

## Rules

- Keep one current guide per topic and link instead of duplicating it.
- Keep the README focused on project purpose, capabilities, installation and a
  first example. Put detailed YAML/backend/API instructions in their user
  references and link to them from the README and documentation index.
- Keep current capability status in the feature matrix and dated execution or
  downstream platform acceptance in reports. User guides retain the scientific
  limits needed to configure a calculation, without repeating platform delivery
  disclaimers or historical run logs.
- User pages describe installation, configuration, running and results.
  Maintainer commands, CI configuration and site operations belong in
  `docs/developer/`.
- Reports are evidence. Link them from the feature matrix or release notes when
  needed; do not require users to read them.
- Keep reports for the current support line in `docs/reports/`; move older
  dated reports to `docs/archive/reports/`.
- Keep one release note per version. Older notes may be archived after the
  current release and its two predecessors.
- Active Markdown links must resolve inside the repository. Use relative links.
- Generated documents must state their source; for example
  `docs/user/YAML_INPUT_VARIABLES.md` is generated from `config_schema.py`.
- Program output and new diagnostics/docstrings use English. Existing
  documentation language is preserved.
- Do not add metadata headers, L1-L4 report levels, a separate report ledger or
  a pending-delete inventory. Use `docs/archive/` directly.

## Testing from a worktree

If the active environment has `atst_tools` installed from another checkout,
run tests from the worktree with either:

```bash
PYTHONPATH=src python -m pytest tests/unit -q
```

or install the worktree in editable mode:

```bash
python -m pip install -e .
```

Verify that `atst_tools.__file__` points to the current worktree before
interpreting test failures.

## Checks

Documentation-only changes should run:

```bash
git diff --check
python scripts/check_docs_governance.py
python -m pytest tests/unit/test_docs_governance.py -q
```

The governance script checks conflict markers, relative links, required entry
documents and basic HTML report parsing. Release facts, report ledgers and
plan registration are no longer separate documentation gates.
