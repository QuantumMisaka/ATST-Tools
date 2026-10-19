# Constant-potential develop-branch smoke evidence (2026-10-20)

Two site builds were exercised on the Pt/LCAO fixture:

| Directory | Engine | Job |
| --- | --- | --- |
| `reference_results.json`, `compensated_results.json` | `develop-git-079fd0c-260724-sm70-auto` (`v3.11.0-beta6`, latest-format reader) | `1864431` |
| `beta9/reference_results.json`, `beta9/compensated_results.json` | `abacus/v3.11.0-beta9-auto` (`v3.11.0-beta9`) | `1864468` |

Both jobs completed with `reference_fcp` and `compensated_gate` result status
`complete`. The `*.log` files are the ATST CLI logs; `slurm-*.out` and
`sacct.txt` record the allocation and completion.
