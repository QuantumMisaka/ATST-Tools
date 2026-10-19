# ABACUSLite Wrapper Guide

## Current Integration State

ATST-Tools uses the independently installed `abacuslite>=1.2.1,<2` package
for every ABACUS calculator and file reader/writer. Missing or broken imports
fail visibly. The old `atst_tools.external.ASE_interface.abacuslite` import paths
forward to the same installed implementation for existing callers, including
Paimon; they do not select a separate backend. 2.3.0 introduced this
required-dependency cutover; 2.3.1 requires `abacuslite>=1.2.1` for the
entropy-aware `free_energy` handoff.

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
checks; 1.2.1 adds the entropy-aware `free_energy` fix, and its baseline did
not repeat the earlier engine validation; upgrading a deployed runtime still
requires its own validation.

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

### Force / fmax 口径（RAW vs 投影，2026-08-04 确认）

ABACUS 的 `running_*.log` 中 `TOTAL-FORCE` 输出的是**全原子原始力（RAW）**，包含
STRU 中 `m 0 0 0` 固定原子——这是期望行为（与 VASP OUTCAR 等 DFT 惯例一致），
固定原子由驱动层约束处理。下游消费口径：

- **计算器层（独立 abacuslite backend）**：返回 RAW 全原子力，不做预投影；约束投影由 ASE
  `Atoms.get_forces()`（默认 `apply_constraint=True`，按 `FixAtoms`）执行。
- **优化器层（sella/ASE）**：收敛判据使用**投影后自由原子力**；不要把 RAW 全原子
  fmax 直接与收敛阈值比较。
- **汇总/报告层（transition summary 等）**：报告 fmax 须明确口径（投影后自由原子力
  与 RAW 全原子力），约束体系收敛以投影后为准。
- **诊断/取证**：比较力必须同口径（投影后 vs 投影后，或仅自由原子原始力）；禁止把
  ASE 投影后值与 running log RAW 全原子值直接对比下结论。

该口径是"固定 Au 上 sella 提前收敛"曾被误诊为力读取 bug 的根因（约束投影被误读为
陈旧缓存），已按此修正。

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
