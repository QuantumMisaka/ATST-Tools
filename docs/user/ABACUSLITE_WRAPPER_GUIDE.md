# ABACUSLite Wrapper Guide

## Current Integration State

ATST-Tools uses the independently installed `abacuslite>=1.2.1,<2` package
for every ABACUS calculator and file reader/writer. Missing or broken imports
fail visibly. The old `atst_tools.external.ASE_interface.abacuslite` import paths
forward to the same installed implementation for existing callers. Version
1.2.1 supplies the explicit `energy_without_entropy` metadata used for
reference-only constant-potential comparisons; see the
[energy definitions](CONFIG_REFERENCE.md#211-constant-potential-evaluation).

This means ABACUS-backed workflows can run through ASE calculators without
requiring users to write Python calculator boilerplate. The supported workflow
entry point remains:

```bash
atst run CONFIG.yaml
```

`neb`, `autoneb`, `d2s`, `dimer`, `sella`, `ccqn`, `relax`, `vibration`, `irc`, and
`md` can all use `calculator.name: abacus`.

### SCF frame identity

For `calculation: scf`, the reader automatically compares the current STRU with
every parsed running-log frame in the same ASE atom order and Cartesian Å
coordinates. It accepts periodic images: a frame is a match when
each atom differs by an integer combination of the same 3×3 cell vectors,
including different integer translations for different atoms and non-orthogonal
cells. The reader scans from the end and returns the latest matching frame, so
an accumulated running log does not require the current structure to be its
last frame.

The check remains fail-closed. Different cells, atom counts or shapes,
non-finite or singular cells/coordinates, element or atom-order mismatches,
malformed frames, and non-integer residual displacements are rejected with a
diagnostic rather than being silently accepted as the current structure. The
matching step does not wrap or modify the STRU, input or parsed `Atoms`, forces,
cell, or transition-state/NEB coordinates.

This automatic behavior has no `pbc_mode`, tolerance, or other YAML switch.
Native `relax`, `md`, and `MD_dump` output keeps its existing last-frame
semantics; the PBC-aware identity check applies only to SCF result selection.
When `write_input` uses a custom `stru_file`, the same filename is used for the
identity check; without a prior `write_input`, the default `STRU` behavior is
preserved. These semantics are part of the current installed abacuslite baseline;
ATST does not patch the parser. The 2.3.0 closeout records installed-package
checks. Version 1.2.1 preserves standard ASE energies and adds separate
no-entropy metadata. Package-level verification is distinct from the
validation required when upgrading a deployed runtime.

## Historical Snapshot Evidence

The former vendored implementation was compared against the independent
package before this cutover. Its frozen baseline and patch ledger remain in
`ABACUSLITE_SNAPSHOT.md` and `PATCHES.md`; they are audit history, not a live
fallback or a current parser contract.
ABACUS parser and calculator changes now belong in the independent
`QuantumMisaka/abacuslite` repository. The old import paths are thin
forwarders, and the ATST consumer suite verifies their object identity.

ATST workflow YAML already defaults NEB, AutoNEB, and the D2S rough DyNEB path
to ASE's recommended `improvedtangent` method. Direct `AbacusNEB(...)`
construction pins the same default explicitly, so its behavior does not
depend on an ASE release's implicit default.

## Forces and fmax: raw and constrained values

ABACUS `TOTAL-FORCE` in `running_*.log` contains **raw forces on all atoms**,
including atoms fixed with `m 0 0 0` in STRU. The driving workflow applies
constraints; the raw log alone does not describe optimizer convergence.

- **Calculator:** the independent abacuslite backend returns raw forces on all
  atoms. ASE `Atoms.get_forces()` applies constraints such as `FixAtoms` by
  default (`apply_constraint=True`).
- **Optimizer (Sella/ASE):** convergence uses forces after constraints are
  applied. For fixed atoms, compare the free-atom fmax with the threshold.
- **Summary/report:** label whether fmax refers to raw all-atom forces or forces
  after constraints. Use the latter to assess constrained optimization.
- **Diagnostics:** compare like with like: constrained forces against constrained
  forces, or raw forces on free atoms only. Comparing an ASE constrained fmax
  directly with the raw all-atom log fmax can misidentify a force-reading error.

This distinction, confirmed on 2026-08-04, explained an earlier report of Sella
converging too early on fixed Au atoms: constraint projection had been mistaken
for stale cached forces.

## Wrapper Boundary

ATST-Tools is a layered wrapper around abacuslite:

- It owns YAML validation, schema defaults, workflow dispatch, restart helpers,
  trajectory naming, and common pre/post-processing.
- It uses abacuslite as the ASE calculator backend for ABACUS calculations.
- It exposes conservative ABACUS input/output helpers for repeated user tasks.
- It does not replace ABACUS, abacuslite, a scheduler, or site-specific
  launchers.

## Common Commands

Inspect the normalized config before launching a run:

```bash
atst config validate config.yaml --print-normalized
```

Prepare ABACUS input files from the ABACUS calculator block:

```bash
atst abacus prepare config.yaml --structure inputs/init.stru --output-dir abacus_input
```

This writes:

- `INPUT` from `calculator.abacus.parameters`.
- `KPT` from `calculator.abacus.kpts` or `parameters.kpts`.
- `STRU` from the supplied structure and `pseudopotentials` / `basissets`.

Collect a conservative output summary:

```bash
atst abacus collect run_abacus --output abacus_results.json
```

The summary records detected `INPUT`, `KPT`, `STRU`, and `running*.log` files.
When the directory contains the files required by the active abacuslite reader,
the command parses the final frame and can export it:

```bash
atst abacus collect run_abacus --output abacus_results.json --structure final.extxyz
```

The collector copies parse inputs into a temporary directory before invoking
abacuslite readers, so original ABACUS outputs are not modified.

## Complex Workflows

Complex workflows are still launched with `atst run`:

```bash
atst run examples/01_neb_Li-Si/config.yaml
atst run examples/08_d2s_Cy-Pt/config.yaml
atst run examples/10_irc_H2/config.yaml
```

For D2S, ATST-Tools uses the same ABACUS calculator backend through rough NEB
and the selected single-ended method. For NEB and AutoNEB, endpoint
single-point governance repairs missing, placeholder, or unmarked (e.g.
uploaded-chain/foreign) endpoint results before ASE NEB construction when
configured with `endpoint_singlepoint: auto`, and recomputes both endpoints under
`always`; `never` preserves user-provided readable endpoint results.

## NEB Image-Level MPI

ATST-Tools assigns each NEB image an independent `Abacus` calculator
directory and uses ASE image-level parallelism for ABACUS-backed NEB and
AutoNEB.

Use an MPI-enabled Python environment compatible with the site launcher and
the installed ABACUS runtime. For ordinary NEB, the outer MPI world size must
equal the number of interior images; for AutoNEB, it must equal
`calculation.n_simul`. Start with `calculator.abacus.mpi: 1`; increasing it
adds a second, inner MPI layer for each ABACUS image calculation.

The outer launcher remains outside ATST-Tools. Configure the ABACUS subprocess
separately with `calculator.abacus.command`. When a bare single-process
`abacus` command runs inside an image-level MPI workflow, ATST-Tools clears the
outer MPI launcher environment for the ABACUS subprocess so it remains a
one-image calculation.

Before launching an image-parallel calculation, check both the import and a
two-rank launch in that environment:

```bash
python -c "from mpi4py import MPI; print(MPI.Get_library_version())"
mpiexec -n 2 python -c "from mpi4py import MPI; print(MPI.COMM_WORLD.Get_rank())"
```

The second command should report ranks `0` and `1`. Serial users do not need
`mpi4py` or these checks.

If the installed wheel is incompatible with the site MPI implementation, rebuild
`mpi4py` with the site compiler:

```bash
MPICC="$(command -v mpicc)" \
  python -m pip install --no-cache-dir --force-reinstall --no-binary=mpi4py mpi4py
```

The site's ABACUS executable and outer launcher must use the same MPI implementation.

For site setup, example validation, and maintainer operations, use the
[developer handover](../developer/HANDOVER.md) and the
[example validation operations guide](../developer/EXAMPLE_VALIDATION_OPERATIONS.md).

## Non-Goals

- No scheduler submission command is provided in this layer.
- No site-specific environment setup is encoded in YAML.
- No complete ABACUS output database is built from run directories.
- No guarantee is made that every abacuslite IO function is exposed at the CLI.

Future expansion should keep this boundary: add helpers only when they reduce a
repeated ATST workflow step and can be tested without launching expensive
calculations.
