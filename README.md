# ATST-Tools

[![Version](https://img.shields.io/badge/version-2.3.0-blue)](pyproject.toml)
[![Python](https://img.shields.io/badge/python-%3E%3D3.10-blue)](pyproject.toml)
[![License](https://img.shields.io/badge/license-LGPL--v3-blue)](#license)

> User documentation and examples:
> [Chinese user guide: GitHub](https://github.com/QuantumMisaka/atst-tools/blob/main/docs/user/USER_GUIDE_CN.md)
> · [Gitee](https://gitee.com/jamesmisaka/atst-tools/blob/main/docs/user/USER_GUIDE_CN.md)
> · [Examples: GitHub](https://github.com/QuantumMisaka/atst-tools/tree/main/examples)
> · [Gitee](https://gitee.com/jamesmisaka/atst-tools/tree/main/examples)
>
> GitHub is the source of truth for development and specifications. Gitee is a
> read-only mirror synchronized manually by the maintainer; PyPI publishes only
> the Python package.

ATST-Tools is a pip-installable **ASE transition-state workflow toolkit** for
ABACUS and DeePMD-kit calculators. It turns the project's legacy script
collection into one governed command-line interface:

```bash
atst run CONFIG.yaml
```

Use it when you want repeatable NEB, AutoNEB, Dimer, Sella, CCQN, D2S,
relaxation, vibration, IRC, MD, or experimental DMF calculations driven by YAML
or by an Agent's Python integration instead of one-off scripts.

## At A Glance

| Area | Current 2.3.0 status |
| :--- | :--- |
| Package | Python package with the `atst` console command. |
| Main interface | `atst run CONFIG.yaml` for all calculator-backed workflows. |
| Lightweight tools | `atst config`, `atst abacus`, `atst neb`, `atst traj`, `atst dimer`, `atst relax`, `atst vibration`. |
| Calculators | ABACUS through `abacuslite`; DeePMD-kit through `deepmd.calculator.DP`. |
| Configuration | Pydantic-governed YAML schema with generated user documentation. |
| Release | `2.3.0` is published on PyPI (2026-10-19). See the [2.3.0 release notes](docs/releases/RELEASE_NOTES_2.3.0.md). |

Version 2.3.0 makes the independently maintained `abacuslite` package the
single ABACUS ASE backend. The previous vendored snapshot is no longer a
runtime fallback.

## What You Can Run

| `calculation.type` | Workflow | Notes |
| :--- | :--- | :--- |
| `neb` | NEB / DyNEB | Endpoint single-point governance is enabled by default. |
| `autoneb` | AutoNEB | Adaptive image insertion plus final-chain post-processing. |
| `d2s` | Double-ended to single-ended TS search | Rough NEB followed by Dimer, Sella, or CCQN. |
| `dimer` | ASE Dimer | Single-ended transition-state search. |
| `sella` | Sella saddle search | Uses the external `sella` package. |
| `ccqn` | CCQN saddle search | Cone-shaped constrained quasi-Newton TS optimization. |
| `relax` | Structure optimization | ASE optimizer based relaxation. |
| `vibration` | Vibrations and thermochemistry | Harmonic and ideal-gas helpers. |
| `irc` | Sella IRC | Sella-backed IRC orchestration with controlled boundary diagnostics. |
| `md` | Molecular dynamics | ASE-driven MD with ABACUS/DP calculators, or ABACUS native MD input/run/output orchestration. |
| `dmf` | Direct MaxFlux | Experimental TS candidate/path optimizer; non-periodic first and not a validated TS result. |

`constant_potential` remains a development candidate and does not establish
production readiness or Paimon/public tool-chain acceptance. Its current scope
is ABACUS-only fixed-geometry single-point/serial scans plus fixed-cell
`relax`/`neb` with the explicit `compensated_gate` boundary. The legacy
`reference_fcp` boundary is reference-only and is rejected for optimization
workflows. See the [configuration reference](docs/user/CONFIG_REFERENCE.md)
before using it.

Local pre/post-processing commands are intentionally lightweight. They do not
construct calculators or submit expensive calculations:

```bash
atst neb make ...
atst neb post ...
atst traj collect ...
atst traj transform ...
atst dimer make-from-neb ...
atst relax post ...
atst vibration post ...
atst config validate ...
atst abacus prepare ...
atst abacus collect ...
```

## Installation

### From PyPI (Recommended)

```bash
pip install atst-tools
```

This installs the current PyPI release. ATST-Tools requires Python 3.10 or
newer; `sella>=2.5` and `abacuslite>=1.2.0,<2` are required runtime
dependencies.

Installation profiles:

- **Serial (default)**: `pip install atst-tools`; installs no `mpi4py` and
  supports every non-image-parallel workflow, including serial NEB/AutoNEB.
- **DP**: `pip install "atst-tools[dp]"`; required only for DeePMD-kit
  calculator workflows.
- **Image-parallel NEB/AutoNEB**: `pip install "atst-tools[parallel]"`; required
  only when launching one Python rank per active image through an external
  `mpirun`/`srun`.

### From Source

```bash
git clone https://github.com/QuantumMisaka/atst-tools.git
cd atst-tools
pip install .
```

Maintainers who need development installs, release artifacts, or validation
operations should use the [developer handover](docs/developer/HANDOVER.md) and
the [example validation operations guide](docs/developer/EXAMPLE_VALIDATION_OPERATIONS.md).

ATST-Tools itself installs the Python workflow layer. Real calculations also
need the selected calculator runtime:

- **ABACUS**: an executable ABACUS installation plus pseudopotential/orbital
  files referenced by YAML.
- **DP / DeePMD-kit**: install `atst-tools[dp]` or provide a compatible
  DeePMD-kit Python installation, plus a model file outside git-tracked paths.
- **MPI image parallelism**: install `atst-tools[parallel]` in an MPI-compatible
  Python environment and launch ATST with one Python rank per active image.
- **DMF**: ATST-Tools vendors PyDMF, but runtime still requires `cyipopt` and
  IPOPT, for example from conda-forge.

## Quick Start

Choose the path that matches what you need:

| Need | Start here |
| :--- | :--- |
| Learn the CLI and inputs in 10 minutes | [Chinese user guide](docs/user/USER_GUIDE_CN.md) |
| Pick an example | [Examples overview](examples/README.md) |
| Check supported features | [Feature status matrix](docs/reports/FEATURE_STATUS_MATRIX.md) |
| Look up YAML semantics | [Configuration reference](docs/user/CONFIG_REFERENCE.md) |
| Look up every schema field | [YAML input variables](docs/user/YAML_INPUT_VARIABLES.md) |
| Use CLI helper commands | [CLI reference](docs/user/CLI_REFERENCE.md) |
| Embed a workflow in Python | [Stable Python API reference](docs/user/PYTHON_API_REFERENCE.md) |
| Browse all documentation paths | [Documentation index](docs/index.md) |

After installing the package, obtain the examples separately if you do not
already have a checkout. Use the tag matching your installed version:

```bash
git clone --branch main https://github.com/QuantumMisaka/atst-tools.git
cd atst-tools
atst config validate examples/06_relax_H2-Au/config.yaml --print-normalized
atst run --dry-run examples/06_relax_H2-Au/config.yaml
```

Try a local operation without ABACUS, a DP model, or a GPU:

```bash
cd examples/09_lightweight_cli
atst neb make inputs/init.xyz inputs/final.xyz 3 -o first_chain.traj --method linear
atst neb summary inputs/neb_result.extxyz --n-max 1 --tail 5
cd ../..
```

`first_chain.traj` contains five initial-guess frames, not evaluated energies or
forces; configuration validation does not establish backend readiness. Once
ABACUS and its data files are available, run the H2/Au Relax example (four
backend MPI ranks, `cusolver`, one-step optimization limit) from its directory
and inspect convergence separately. The [Chinese guide](docs/user/USER_GUIDE_CN.md)
walks through that run, its outputs, and common problems.

Embed the same schema validation from Python when a calling program needs a
structured result rather than terminal output:

```python
from atst_tools.api import validate_config

config = validate_config("examples/06_relax_H2-Au/config.yaml")
print(config["calculation"]["type"])
```

CLI/YAML and the Agent Python API use the same ASE transition-state workflows.
Use CLI/YAML for interactive or scheduled calculations; use the API for
structured results or caller-owned calculator integration. A calculation does
not by itself establish transition-state validity or launch IRC. Validate a
candidate separately with `atst vibration validate` or the additive
`validate_ts` API, which consume saved vibration and optional frozen-structure
evidence. The [Python API reference](docs/user/PYTHON_API_REFERENCE.md) defines
the result, MPI, artifact, and calculator boundaries.

For an external process host that needs a stable JSON handoff without parsing
terminal output, use the installed API runner:

```bash
python -m atst_tools.api.runner --config config.yaml --workdir run --result-json atst_api_result.json
```

It invokes the same configuration-driven API and writes `atst-api-result-v1`
only from root rank. ATST does not allocate resources or start the outer Python
MPI launch; a configured calculator command may launch its own ABACUS MPI process.
See the [Python API
reference](docs/user/PYTHON_API_REFERENCE.md#process-runner-for-external-hosts)
for flags, exit codes, and artifact ownership.

Print a schema-governed template:

```bash
atst run --show-template neb --calculator abacus
atst run --show-template neb --calculator dp
atst run --show-template dmf --calculator dp
```

List available workflow types:

```bash
atst run --list-types
```

Print the project banner and contributor credits:

```bash
atst banner
```

## Minimal YAML Shape

Every production workflow uses the same top-level structure:

```yaml
calculation:
  type: neb
  init_chain: inputs/init_neb_chain.traj
  fmax: 0.05
  max_steps: 100
  climb: true
  two_stage: true
  stage1_steps: 20
  stage1_fmax: 0.20

calculator:
  name: abacus
  abacus:
    command: abacus
    mpi: 4
    omp: 1
    directory: run_neb
    kpts: [2, 2, 2]
    parameters:
      calculation: scf
      ecutwfc: 100
      basis_type: lcao
```

The governed schema defines defaults, types, and descriptions for user-facing
inputs. See [YAML input variables](docs/user/YAML_INPUT_VARIABLES.md) for the
generated reference and [configuration reference](docs/user/CONFIG_REFERENCE.md)
for hand-written guidance.

## Calculator Backends

### ABACUS

ABACUS is integrated through the independently installed `abacuslite` package.
In 2.3.0 it is a required dependency; missing or broken imports
fail visibly. Calculator, INPUT/STRU/KPT helpers, output readers, reverse
configuration, native MD and gate compensation use that one package. Historical
imports under `atst_tools.external.ASE_interface.abacuslite` remain as forwarding
paths for existing consumers and do not contain a second implementation.

Typical example files use:

```yaml
calculator:
  name: abacus
  abacus:
    command: abacus
    mpi: 4
    omp: 1
    pseudo_dir: ../data
    orbital_dir: ../data
```

ATST-Tools is a layered `abacuslite` wrapper. Calculator-backed workflows such
as NEB, D2S, Dimer, Sella, CCQN, Relax, Vibration, and IRC still run through
`atst run CONFIG.yaml`; local ABACUS helpers support safe input preparation and
result collection:

```bash
atst abacus prepare config.yaml --structure inputs/init.stru --output-dir abacus_input
atst abacus collect run_neb --output abacus_results.json
```

These helper commands do not run ABACUS and do not submit scheduler work.

### DeePMD-kit / DP

DP support uses the official ASE calculator entry point,
`deepmd.calculator.DP`. DeePMD-kit detects the model backend from the model
file. Multi-head DPA/DPA3 models are configured through `calculator.dp.head`.

```yaml
calculator:
  name: dp
  dp:
    model: ../../temp_repos/dp_model/DPA-3.1-3M.pt
    head: Omat24
    omp: 4
    share_calculator: true
```

The 2.0.0 DP validation used DPA-3.1-3M with the `Omat24` head. The pinned
download source, checksum, expected size, and local path are recorded in
`examples/dp_model_manifest.json`. Model files and runtime outputs are
intentionally not tracked by git.

```bash
python scripts/download_dp_model.py
python scripts/download_dp_model.py --check-only
```

## Examples

The [examples directory](examples/README.md) is the fastest way to learn the
project. Each calculation example uses `config.yaml` for ABACUS and, where
available, `config_dp.yaml` for DP.

- **NEB / AutoNEB / D2S**: `01_neb_Li-Si`, `02_neb_H2-Au`,
  `03_autoneb_Cy-Pt`, `08_d2s_Cy-Pt`, `13_neb_parallel_Cy-Pt`,
  `14_autoneb_parallel_Cy-Pt`.
- **Single-ended TS search**: `04_dimer_CO-Pt`, `05_sella_H2-Au`,
  `12_ccqn_H2-Au`.
- **Relaxation, vibration, IRC**: `06_relax_H2-Au`, `07_vibration_H2-Au`,
  `10_irc_H2`, `11_vibration_ideal_gas_H2`.
- **MD and local CLI**: `15_md_Li-Si`, `09_lightweight_cli`.
- **Constant potential**: `19_constant_potential_Pt` (development fixtures, not
  a production benchmark).

## For Developers

Start with the [source architecture guide](src/atst_tools/README.md) to understand
package responsibilities, the shared CLI/API execution path, the purpose of
`runtime` and `bench`, and where to make changes.

Developer governance starts from these maintained entry points:

- [Documentation index](docs/index.md)
- [Handover checklist](docs/developer/HANDOVER.md)
- [YAML input governance](docs/developer/YAML_INPUT_GOVERNANCE.md)
- [Documentation guide](docs/developer/DOCUMENTATION.md)
- [ABACUS wrapper guide](docs/user/ABACUSLITE_WRAPPER_GUIDE.md)
- [Maintained atst-cli skill](docs/skills/atst-cli/SKILL.md)
- [Release readiness gate](docs/developer/PYPI_RELEASE_AUTOMATION.md#local-readiness-gate)

Project status and validation entry points:

- [Feature status matrix](docs/reports/FEATURE_STATUS_MATRIX.md)
- [2.3.0 release notes](docs/releases/RELEASE_NOTES_2.3.0.md)
- Historical validation reports are indexed in
  [the archive](docs/archive/README.md).

## Version Governance

The package version has one source of truth: `pyproject.toml ->
[project].version`. Runtime entry points read it through
`atst_tools.package_version()`; source-tree runs read `pyproject.toml`, while
installed runs use distribution metadata generated from the same field. User
YAML has no `config_version` and is governed by the installed package schema;
unknown top-level fields are rejected.

## Project Boundary

ATST-Tools owns workflow orchestration, YAML validation, calculator construction,
trajectory naming, restart handling, ABACUS input/output helpers, examples, and
documentation. Numerical engines remain external:

- ABACUS owns first-principles electronic-structure calculations.
- DeePMD-kit owns DP model loading and inference.
- ASE owns the core optimizer and transition-state method implementations.
- Sella owns its saddle-search and IRC algorithms.

## References

If you publish work that uses the transition-state workflows in ATST-Tools,
please cite the underlying methods alongside ATST-Tools itself:

- ATST-Tools and the double-to-single (D2S) transition-state workflow:
  Liu, Z.-Q.; Deng, Z.; Zhao, H.; Wang, H.; Chen, M.; Jiang, H.
  *FT²DP: Large Atomic Model Fine-Tuned Machine Learning Potential for
  Accelerating Atomistic Simulation of Iron-Based Fischer–Tropsch Synthesis.*
  J. Mater. Inform. **5** (2), 27 (2025).
  <https://doi.org/10.20517/jmi.2024.105>
- Sella saddle-point search:
  Hermes, E. D.; Sargsyan, K.; Najm, H. N.; Zádor, J.
  *Sella, an Open-Source Automation-Friendly Molecular Saddle Point Optimizer.*
  J. Chem. Theory Comput. **18** (11), 6974-6988 (2022).
  <https://doi.org/10.1021/acs.jctc.2c00395>
- CCQN (cone-shaped constrained quasi-Newton) transition-state optimization:
  Wu, Y.; Wang, H. *Cone-Shaped Constrained Quasi-Newton Method: Efficient
  and Robust Single-Ended Transition State Optimization Algorithm.*
  J. Chem. Theory Comput. **21** (18), 9054-9065 (2025).
  <https://doi.org/10.1021/acs.jctc.5c01015>

## License

LGPL-v3 License.
