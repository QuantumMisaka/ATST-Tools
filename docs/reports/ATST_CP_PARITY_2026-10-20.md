# ATST constant-potential reference parity (2026-10-20)

**Status**: Historical candidate single-point `reference_fcp` parity passed
on a real ABACUS engine. The release handoff was corrected during review;
see the closeout note below for its distinct verification scope.

## Scope

This is an integration-reproducibility check, not a scientific convergence or
production-accuracy acceptance. It compares:

- FCP v2 (the external `FCP-abacus-ase` reference implementation) with
  `ase-abacus` and a real ABACUS LTS executable.
- ATST `2.3.0` source with the entropy-aware energy mapping and `abacuslite`.

The case is the Pt LCAO slab from
`examples/19_constant_potential_Pt`, materialized with 217 initial electrons,
`N0=216`, `U=0.4 V`, `work_ref=4.6 eV`, Gamma `1x1x1`, `ecutwfc=100 eV`,
`scf_thr=1e-9`, and `cal_stress=1`.

## Versions and environment

| Item | Value |
| --- | --- |
| SAI job | `1864305`, `4V100`, 1 GPU, `COMPLETED`, 00:03:22 |
| ABACUS | LTS `3.10.1-sm70-auto` |
| Python / ASE | 3.13.15 / ASE 3.29.0 for ATST; ase-abacus `3.26.0b1` for FCP v2 |
| FCP v2 source | `version2/FCPelectrochem.py`, SHA256 `bd1cbf41abf34d68a41f6f1579275e6c082e8f0f8be09c2d629153c2aa745040` |
| ATST | current source (2.3.0 line plus CP energy-mapping fix) |
| abacuslite | 1.2.0 source plus the `free_energy` fix; local wheel SHA256 `ad05df5afd63be6298a6392a40a077002f3284f32af84435336b1602d31abd67` |
| Tolerances | `FCPconv=1e-4`, `potential_tolerance_v=1e-4`, `always_adjust=false`, `adaptive_lr=false` |

The published `abacuslite==1.2.0` does **not** contain this local fix. The
parity result applies to the patched wheel used in this run.

## Procedure

1. `scripts/run_fcp_v2_reference.py` runs FCP v2 and writes
   `fcp_reference.json`.
2. `atst_tools.scripts.cli run config_reference.yaml` runs ATST with
   `energy_boundary: reference_fcp`.
3. `scripts/extract_constant_potential_record.py` converts the ATST point into
   the same record schema.
4. `scripts/compare_constant_potential_results.py` compares the records.

The committed evidence is in
`docs/reports/data/ATST_CP_PARITY_20261020/`.

## Result

`comparison.json` reports `ok: true` with `mismatches: []`.

| Quantity | FCP v2 | ATST | Absolute delta |
| --- | ---: | ---: | ---: |
| `nelec` | 215.966013 | 215.966013788 | 7.9e-7 |
| `efermi_eV` | 5.034403 | 5.034404288 | 1.3e-6 |
| `vacuum_level_eV` | 10.034402 | 10.034401012 | 9.9e-7 |
| `mu_calc_eV` | -5.0 | -4.999996724 | 3.3e-6 |
| `residual_mu_eV` | 0.0 | 3.2762e-6 | 3.3e-6 |
| `energy_eV` | -39571.255745341 | -39571.255745094 | 2.5e-7 |
| `free_energy_eV` | -39572.373860757 | -39572.373860605 | 1.5e-7 |
| maximum force component | — | — | 2.7e-7 eV/A |

Informational differences (not part of the pass/fail decision):

- Newton iteration count: FCP v2 `6`, ATST `5`.
- Capacitance estimate: FCP v2 `0.199198052`, ATST `0.198502026`
  (delta `6.96e-4 e/V`). The two implementations use different fits.

## Finding and fix

The first real-engine attempt exposed a 1.118 eV difference in the CP
`energy` field. The ABACUS log showed:

```text
E_KohnSham     -39572.545121
E_entropy(-TS)  -1.118110
```

FCP v2 publishes the no-entropy energy (`Ewithoutentropy`) under its `energy`
key and `E_KohnSham` under `free_energy`. Published abacuslite 1.2.0 reports
both standard keys as `E_KohnSham`, so it lacks the separate value needed to
reproduce that reference bookkeeping.

The initial parity candidate put the no-entropy value in backend `free_energy`
and swapped the keys in ATST. **That backend mapping was rejected during
pre-release review:** ABACUS `E_KohnSham` already includes `E_entropy(-TS)`;
subtracting it is not the ASE force-consistent free energy. The unconditional
swap also affected the compensated-gate optimization path.

The corrected 1.2.1 / 2.3.1 pair preserves backend
`energy=free_energy=E_KohnSham`, exposes `energy_without_entropy` separately,
and uses that extra value only for `reference_fcp.energy`.
`compensated_gate` retains the standard backend keys and correction.

### Release closeout verification (2026-10-20)

The retained SAI outputs below are unchanged historical evidence for their
recorded patched wheel. Parser and file-calculator regressions check the
corrected handoff on legacy/latest fixtures and matching-frame SCF selection.
ATST regressions distinguish both CP boundaries. The test
`test_recorded_lts_trace_keeps_reference_parity_with_explicit_energy` replays
all five retained electronic evaluations through the corrected CP calculator
and checks the final two energy values against FCP v2 within `1e-6 eV`.
This is offline replay, not a new GPU run; the develop compensated-gate values
below must not be represented as corrected-release results.

## What this confirms

- The ATST `reference_fcp` single-point path is numerically consistent with
  FCP v2 for the tested Pt/LCAO case on ABACUS LTS 3.10.1.
- The abacuslite handoff (`efermi`, vacuum/fermishift, forces, and now
  entropy-aware energy) is sufficient for this case.

## Develop/latest-format smoke

Two site develop builds were also exercised on the same Pt/LCAO fixture:

| Job | Engine | `reference_fcp` | `compensated_gate` |
| --- | --- | --- | --- |
| `1864431` | `develop-git-079fd0c-260724-sm70-auto` (`v3.11.0-beta6`) | complete | complete |
| `1864468` | `abacus/v3.11.0-beta9-auto` (`v3.11.0-beta9`) | complete | complete |
| `1864773` | locally built `origin/develop` `7be5df975` (`v3.11.0-beta11`) | complete | complete |

Key values:

| Case | `nelec` | `efermi` (eV) | `mu_calc` (eV) | residual (eV) | `energy` (eV) |
| --- | ---: | ---: | ---: | ---: | ---: |
| beta6 `reference_fcp` | 215.966013787 | 5.034404286 | -4.999996729 | 3.27e-6 | -39571.255745095 |
| beta9 `reference_fcp` | 215.966013788 | 5.034404287 | -4.999996725 | 3.28e-6 | -39571.255745099 |
| beta6 `compensated_gate` | 217.0 | 9.325994439 | 15.045970658 | -1.77e-10 | -39581.254151944 |
| beta9 `compensated_gate` | 217.0 | 9.325994439 | 15.045970658 | 7.40e-11 | -39581.254151936 |
| beta11 `reference_fcp` | 215.966013526 | 5.034404363 | -4.999997234 | 2.77e-6 | -39571.255883550 |
| beta11 `compensated_gate` | 217.0 | 9.325994744 | 15.045970064 | -5.94e-7 | -39581.254287417 |

`abacuslite.core._uses_legacy_io("v3.11.0-beta6")` returns `False`, so the
develop runs use `latestio`. The develop smoke is an ATST-only execution and
artifact check, not an FCP v2 numerical parity run. The beta11 build used
`deepmodeling/abacus-develop` `origin/develop` `7be5df975` with
`ENABLE_LCAO=ON`, `ENABLE_ELPA=ON`, `USE_CUDA=ON`, and
`CMAKE_CUDA_ARCHITECTURES=70`; the binary SHA256 is
`422636c79be9476958045610202fe4414db635b9b3ca8ca9113ba899231859c6`.

## What this does not confirm

- Scientific convergence, k-mesh/cutoff/vacuum convergence, or production
  accuracy.
- `compensated_gate`; FCP v2 is not its oracle.
- FCP v2 numerical parity on the develop build; the LTS parity and develop
  smoke are separate checks.
- Scan, restart, relax, NEB, `always_adjust=true`, or other basis/spin
  combinations.
- The published `abacuslite==1.2.0`, or GPU execution of the corrected
  release pair; see the release notes for indexed-package acceptance.
