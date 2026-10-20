# YAML Input Governance

**Version**: 2.0.0
**Date**: 2026-10-20
**Status**: Maintained
**Owner**: ATST-Tools maintainers

This guide defines how to develop YAML input variables for `atst run`. Every
user-editable variable needs one type, default, meaning, documentation source,
and corresponding test coverage.

## 1. Principles

- **Schema first**: `src/atst_tools/utils/config_schema.py` is the source of
  truth for YAML variables.
- **One canonical path**: Give each user variable one public YAML path. Do not
  introduce aliases, duplicate fields, or compatibility bypasses.
- **Defaults in the schema**: Pydantic must manage defaults except where a
  value is required from the user and has no safe default.
- **Consume normalized configuration**: `atst run` must call
  `ConfigLoader.normalize()` before a workflow or calculator reads fields.
- **Generate the variable table**: Generate the non-calculator YAML variable
  table from the schema; do not maintain it by hand.

User YAML has no `config_version`: the installed package schema validates it,
and unknown top-level fields are rejected.

ABACUS INPUT keywords are the exception. They pass through
`calculator.abacus.parameters` to ABACUS and are not individually defined in
the ATST-Tools schema.

## 2. Source files

- Schema: `src/atst_tools/utils/config_schema.py`
- Loader: `src/atst_tools/utils/config.py`
- Markdown exporter: `src/atst_tools/utils/config_docs.py`
- Generated user variable table: `docs/user/YAML_INPUT_VARIABLES.md`
- Manual user reference: `docs/user/CONFIG_REFERENCE.md`
- Governance tests: `tests/unit/test_config.py`, `tests/unit/test_config_governance.py`

## 3. Add or change a YAML variable

1. Add the field to `config_schema.py` with its type, default, and
   `Field(description=...)`.
2. Add Pydantic validators for enums, positive ranges, mutually exclusive
   fields, and structural constraints as needed.
3. If the field changes calculation behavior, make the workflow or calculator
   consume the normalized value.
4. For internal tests that construct workflows directly, reuse schema defaults
   through `apply_calculation_defaults()`.
5. Regenerate the variable table:
   ```bash
   python -m atst_tools.utils.config_docs --output docs/user/YAML_INPUT_VARIABLES.md
   ```
6. Update `CONFIG_REFERENCE.md`, the relevant user guide, and example YAML.
7. Update unit tests to cover validation/normalization, documentation sync,
   and rejection of redundant variables.

## 4. Disallowed changes

- Do not add user-variable defaults inside a workflow, such as
  `calc_config.get("new_key", 123)`.
- Do not give one semantic variable multiple YAML names.
- Do not place calculator backend-specific parameters in `calculation`.
- Do not edit `YAML_INPUT_VARIABLES.md` by hand without updating the schema
  and exporter.

Dynamic runtime defaults remain valid. For example, MPI `world.size` can
determine `n_simul` when the input sets `n_simul: null`.

## 5. Add an `atst run` feature

The minimum integration path for a new calculation workflow is:

1. Add the implementation under `src/atst_tools/workflows/` or
   `src/atst_tools/mep/`.
2. Add a new `Calculation` schema to `config_schema.py`,
   `VALID_CALCULATION_TYPES`, and the discriminated union.
3. Add `atst run` dispatch and a `--show-template` template in
   `src/atst_tools/scripts/main.py`.
4. Add example YAML and user documentation.
5. Add workflow, configuration, and example unit tests.
6. Regenerate `docs/user/YAML_INPUT_VARIABLES.md`.

For a new calculator backend, integrate through
`src/atst_tools/calculators/factory.py` and put backend-specific fields under
`calculator.<backend>`.
