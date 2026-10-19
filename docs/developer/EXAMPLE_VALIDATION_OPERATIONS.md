# Example Validation Operations

**版本**: 2026-07-23
**日期**: 2026-07-23
**状态**: 维护
**责任人**: ATST-Tools maintainers

## Scope

This guide is the sole active developer entry for maintaining repository
examples, curated outputs, and execution evidence. It is not a user quick
start: users should begin with the runnable learning paths in
[examples/README.md](../../examples/README.md).

Keep the example README limited to ordinary prerequisites, configuration
validation, workflow commands, and scientific limitations. Site setup,
submission patterns, output provenance, and release-facing evidence belong
here. The SAI partition names, QOS names, module names, and job IDs below are
historical evidence for the recorded runs, not public requirements. Reconfirm
current site policy and available resources before submitting new work.

## Local maintainer checks

Start by checking that a changed configuration parses and that its selected
workflow can be planned without launching a calculator:

```bash
atst config validate examples/06_relax_H2-Au/config.yaml --print-normalized
atst run --dry-run examples/06_relax_H2-Au/config.yaml
```

The short DP configurations used for P0/P1 interface coverage are:

```bash
cd examples/02_neb_H2-Au
atst run config_two_stage_dp.yaml

cd ../10_irc_H2
atst run config_descent_dp.yaml

cd ../12_ccqn_H2-Au
atst run config_auto_modes_dp.yaml
```

Their ABACUS counterparts are `config_two_stage.yaml`,
`config_descent.yaml`, and `config_auto_modes.yaml`. Curated outputs for these
cases live under the respective `outputs/` directories. For ordinary NEB, the
main 01/02 examples and example 13 use `two_stage: true` with
`stage1_fmax: 0.20` and `stage1_steps: 20`; the short configurations retain
their recorded step limits.

### Constant-potential FCP v2 reference trace

The committed `tests/data/constant_potential_reference.json` is generated from
the pinned upstream `version2/FCPelectrochem.py` source with a deterministic
analytic backend. Reproduce it without executing ABACUS:

```bash
ATST_FCP_REFERENCE_SOURCE=/path/to/v2_FCPelectrochem.py \
  pytest tests/integration/test_constant_potential_reference_source.py -q
```

The reference source is GPL-3.0 and is not vendored; keep it in a local
reference tree. This checks the FCP v2 outer loop and the committed trace
provenance, not a real-engine calculation. A real-engine side-by-side run
between FCP v2 + ase-abacus and ATST + abacuslite remains an opt-in site
acceptance.

### Constant-potential real-engine record parity

For a real-engine single-point comparison, produce one normalized record from
each implementation and compare them:

```bash
# Candidate side: an ATST run directory containing constant_potential_results.json
python scripts/extract_constant_potential_record.py \
  --input RUN_DIR/constant_potential_results.json \
  --output atst_record.json

# Reference side: requires legacy ase-abacus, pandas, ABACUS, and the pinned
# version2/FCPelectrochem.py checkout
python scripts/run_fcp_v2_reference.py \
  --config fcp_reference_config.json \
  --fcp-source /path/to/v2_FCPelectrochem.py \
  --output fcp_record.json

python scripts/compare_constant_potential_results.py \
  --reference fcp_record.json \
  --candidate atst_record.json
```

The reference producer accepts a JSON configuration with the same physical
settings on both sides. A minimal shape is:

```json
{
  "atoms": "/path/to/STRU",
  "atoms_format": "abacus",
  "command": "mpirun -np 4 abacus",
  "pseudo_dir": "/path/to/PP_ORB",
  "basis_dir": "/path/to/PP_ORB",
  "pseudopotentials": {"Pt": "Pt_ONCV_PBE-1.0.upf"},
  "basissets": {"Pt": "Pt_gga_7au_100Ry_4s2p2d1f.orb"},
  "kpts": [1, 1, 1],
  "parameters": {
    "calculation": "scf",
    "basis_type": "lcao",
    "ecutwfc": 100
  },
  "directory": "/path/to/fcp_reference_run",
  "fcp": {
    "U": 0.4,
    "NELECT": 217.0,
    "NELECT0": 216.0,
    "C": 0.0125,
    "FCPconv": 1e-5,
    "work_ref": 4.6,
    "max_FCP_iter": 100,
    "always_adjust": false,
    "adaptive_lr": false
  }
}
```

Use `always_adjust: false` for the first single-point parity. FCP v2 defaults
to `true`; ATST currently stops at convergence, so the adjusted next guess is
a separate warm-start decision.

The same comparison is available as an opt-in test:

```bash
ATST_FCP_REFERENCE_CONFIG=fcp_reference_config.json \
ATST_FCP_REFERENCE_SOURCE=/path/to/v2_FCPelectrochem.py \
ATST_CP_CANDIDATE_RESULTS=RUN_DIR/constant_potential_results.json \
pytest --run-slurm tests/integration/test_constant_potential_record_parity.py -q
```

The comparator uses explicit absolute tolerances and reports every differing
field. It compares one point; scan, restart, relax/NEB, and cross-basis/spin
matrices remain separate acceptance cases.

The 2026-10-20 Pt/LCAO `reference_fcp` single-point run (SAI job `1864305`,
ABACUS LTS 3.10.1) is recorded in
[ATST_CP_PARITY_2026-10-20.md](../reports/ATST_CP_PARITY_2026-10-20.md). It
passed after the abacuslite entropy-aware `free_energy` fix; the published
`abacuslite==1.2.0` does not yet include that fix.

Develop/latest-format smoke runs on the same fixture are recorded in
`docs/reports/data/ATST_CP_DEVELOP_20261020/`:

- `abacus/develop-git-079fd0c-260724-sm70-auto` (`v3.11.0-beta6`), job
  `1864431`;
- `abacus/v3.11.0-beta9-auto` (`v3.11.0-beta9`), job `1864468`.
- locally built `deepmodeling/abacus-develop` `origin/develop` `7be5df975`
  (`v3.11.0-beta11`), job `1864773`; build and smoke evidence is under
  `docs/reports/data/ATST_CP_DEVELOP_20261020/latest/`.

Both `reference_fcp` and `compensated_gate` completed. These are execution and
format checks, not FCP v2 numerical parity or scientific acceptance.

Use the normal unit and example checks appropriate to the changed surface.
The authoritative P0/P1 execution evidence is
[P0_P1_EXAMPLE_RUNTIME_VALIDATION_2026-05-28.md](../archive/reports/P0_P1_EXAMPLE_RUNTIME_VALIDATION_2026-05-28.md).

## SAI and site-specific execution

Image-level NEB parallelism requires an MPI-aware Python environment. The
recorded SAI pattern loaded `abacus/LTSv3.10.1-sm70-auto`, activated
`atst-neb-mpi`, and used one outer Python rank per active image:

```bash
module load abacus/LTSv3.10.1-sm70-auto
conda activate atst-neb-mpi
cd examples/01_neb_Li-Si
mpirun -np 3 atst run config_parallel_smoke.yaml

cd ../03_autoneb_Cy-Pt
mpirun -np 4 atst run config_parallel_smoke.yaml
```

The outer `mpirun` schedules ASE images. `calculator.abacus.mpi` controls the
ABACUS subprocess count per image; keep it at `1` for an initial parallel
check unless the allocation is sized for nested MPI. Place the outer launch in
the site's submission script, or use the launcher that provides an `mpi4py`
world. Design constraints and validation evidence are in
[MPI4PY_ASE_NEB_PARALLEL_ATST_SUMMARY_2026-05-27.md](../archive/reports/MPI4PY_ASE_NEB_PARALLEL_ATST_SUMMARY_2026-05-27.md).

The dedicated Cy-Pt image-parallel historical submission patterns were:

```bash
cd examples/13_neb_parallel_Cy-Pt
sbatch submit_huge_gpu.sbatch

cd ../14_autoneb_parallel_Cy-Pt
sbatch submit_rush_gpu.sbatch
```

Example 13 recorded five interior images using the `8V100V0` partition and
`huge-gpu` QOS, with outer `mpirun -np 5`, `calculator.abacus.mpi: 4`, and an
isolated inner launcher. Example 14 recorded `n_simul: 4` using `rush-gpu`,
outer `mpirun -np 4`, and `calculator.abacus.mpi: 1`. Both use
`calculator.abacus.omp: 8`. Their input thresholds are respectively
`stage1_fmax: 0.20`, `stage1_steps: 20`, `fmax: 0.12`; and
`fmax: [0.20, 0.20]`.

For experimental DMF staging, run the recorded P3 submission entry only after
the local configuration checks and current site-resource review:

```bash
cd examples/18_dmf_production_validation
export ATST_DMF_ATST=/path/to/atst
sbatch submit_dmf_dp_rush_1gpu.sbatch
```

That script runs two staged DP-backed cases and writes
`dmf_validation_report.json` with `scripts/validate_dmf_candidates.py`. The
same directory contains the recorded D2S-DMF Sella, vibration, descent-IRC,
ABACUS comparison, and refined ABACUS comparison submissions:

- `submit_d2s_dmf_sella_dp_rush_1gpu.sbatch`
- `submit_d2s_dmf_sella_vib_dp_rush_1gpu.sbatch`
- `submit_d2s_dmf_irc_dp_rush_1gpu.sbatch`
- `submit_d2s_dmf_abacus_compare_rush_gpu.sbatch`
- `submit_d2s_dmf_abacus_refined_compare_rush_gpu.sbatch`

DMF-D2S remains experimental: `rough_method: dmf` is opt-in and has not
replaced the default NEB rough stage. When it is selected, ATST uses the final
DMF `t_eval` grid to choose the `tmax` neighbourhood for displacement vectors,
CCQN references, and automatic vibration indices.

## Curated-output provenance

Curated reference values and structures under `examples/` are reviewable
artifacts, not generic user baselines. The ABACUS records are from SAI `4V100`
using ABACUS LTS 3.10.1 with GPU `ks_solver: cusolver`. The DP records are from
SAI `4V100PX` using `temp_repos/dp_model/DPA-3.1-3M.pt` with SHA-256
`86dd3a804d78ca5d203ebf98747e8f16dff9713ba8950097ceb760b161e19907` and head
`Omat24`. The model URL and identity are pinned in
`examples/dp_model_manifest.json`; do not add the checkpoint as a normal Git
blob.

The completed image-parallel outputs have the following historical provenance:

- `13_neb_parallel_Cy-Pt/outputs/`: job `461967` on `8V100V0`; inspect with
  `atst neb summary outputs/neb_parallel_nested_mpi.traj --n-max 5 --strict`.
- `14_autoneb_parallel_Cy-Pt/outputs/`: job `462244` on `4V100`; inspect with
  `atst neb summary --autoneb-prefix outputs/run_autoneb_parallel_single_gpu`.

The retained work directories for the latter runs are under
`validation_runs/cypt_parallel_e2e_20260528/`. Do not overwrite curated inputs
or outputs during ad hoc runs; write generated outputs to a disposable or
ignored directory, then deliberately curate any accepted replacement.

Relevant active evidence reports are:

- [EXAMPLES_MAIN_BRANCH_COMPARISON_LTS3101_2026-05-19.md](../archive/reports/EXAMPLES_MAIN_BRANCH_COMPARISON_LTS3101_2026-05-19.md)
- [DPA3_DP_EXAMPLES_VALIDATION_2026-05-28.md](../archive/reports/DPA3_DP_EXAMPLES_VALIDATION_2026-05-28.md)
- [NEB_IMAGE_PARALLEL_E2E_VALIDATION_2026-05-29.md](../archive/reports/NEB_IMAGE_PARALLEL_E2E_VALIDATION_2026-05-29.md)
- [TWO_STAGE_NEB_LTS3101_VALIDATION_2026-06-04.md](../archive/reports/TWO_STAGE_NEB_LTS3101_VALIDATION_2026-06-04.md)
- [DMF_P3_VALIDATION_STAGING_2026-06-18.md](../archive/reports/DMF_P3_VALIDATION_STAGING_2026-06-18.md),
  [DMF_P3_RUNTIME_VALIDATION_2026-06-18.md](../archive/reports/DMF_P3_RUNTIME_VALIDATION_2026-06-18.md), and
  [DMF_P4_D2S_ABACUS_COMPARISON_2026-06-18.md](../archive/reports/DMF_P4_D2S_ABACUS_COMPARISON_2026-06-18.md)

## Evidence reports and lifecycle

Treat the records in this guide as site-specific historical evidence. When an
example, curated output, or execution boundary changes, update the affected
reference data, write or revise the appropriate evidence report, and update
the feature matrix when the supported status changes. Follow archival rules in
[DOCUMENTATION.md](DOCUMENTATION.md); do not replace a current conclusion by
silently editing an old historical report.

For the complete example and release checklist, return to
[HANDOVER.md](HANDOVER.md). Keep site commands in this guide or evidence
reports, never in the example learning map.
