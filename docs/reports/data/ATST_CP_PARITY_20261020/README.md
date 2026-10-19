# Constant-potential parity evidence (2026-10-20)

SAI job `1864305` (`4V100`, 1 GPU, COMPLETED, 00:03:22), ABACUS LTS
`3.10.1-sm70-auto`, Pt/LCAO `reference_fcp` single point.

| File | Contents |
| --- | --- |
| `fcp_reference.json` | FCP v2 reference record. |
| `atst_candidate.json` | ATST record extracted from the CP result. |
| `comparison.json` | Comparator report; `ok: true`, no hard mismatches. |
| `atst_results.json` | Raw ATST constant-potential result document. |
| `fcp_tmp_log.txt` | FCP v2 `tmp-log-FCP.txt` iteration table. |
| `config_reference.yaml` | ATST input used for the run. |
| `fcp_config.json` | FCP v2 producer input. |
| `slurm-1864305.out` | Slurm job output and final status line. |
| `sacct-1864305.txt` | Slurm accounting record. |

The run used a locally patched `abacuslite` 1.2.0 wheel
(`ad05df5afd63be6298a6392a40a077002f3284f32af84435336b1602d31abd67`)
that exposes the entropy-corrected `free_energy`. The published
`abacuslite==1.2.0` does not include that fix.
