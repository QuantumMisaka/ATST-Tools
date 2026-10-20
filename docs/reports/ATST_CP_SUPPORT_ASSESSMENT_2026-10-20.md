# Constant-potential support assessment — 2026-10-20

The 2.3.2 software supports ABACUS constant-potential single points and serial
scans, plus fixed-cell relaxation, serial NEB and ASE MD under the explicit
`compensated_gate` boundary. This assessment separates that software contract
from numerical validation coverage and downstream deployment.

Ruling: use **Supported** for those implemented workflow contracts and describe
validation coverage separately. The former blanket development-candidate label
hid completed optimization/restart acceptance and conflated application-specific
scientific convergence with software support. This is a documentation reassessment:
no solver, numerical default, supported configuration or persisted-data field
changes. `reference_fcp` retains its reference-only purpose.

## Evidence by capability

| Capability | Software support | Evidence and its scope |
| --- | --- | --- |
| Electronic solve and completed-state handoff | Supported | Regression coverage for measured residuals, SCF failure, electron updates, energy/force conventions, constraints and stale-state rejection. Real LCAO/PW compensated-boundary checks probe energy derivatives and forces. |
| Fixed-geometry single point and serial target scan | Supported | Real GPU single-point/three-target scan with runtime evidence; completed-point persistence and restart regressions. Scan-derived capacitance/PZC conclusions still require adequate sampling. |
| Fixed-cell relaxation | Supported with `compensated_gate` | Historical H/Pt endpoint relaxation and real interruption/continuation met the original force and potential criteria; current restart and constraint regressions retain the software contract. |
| Fixed-cell serial NEB | Supported with `compensated_gate` | Historical constrained H/Pt path and per-image interruption/continuation passed; current endpoint-identity, per-image warm-start and constraint regressions. CP-specific real image-parallel execution remains unverified. |
| Fixed-cell ASE MD | Supported since 2.3.2 with `compensated_gate` | Analytical production-path tests cover all four integrators, force consistency, timestep refinement, constraints, failure and restart. Final reviewed source passed short V100 NVT/Verlet and split-restart runs. Longer real-backend integration and sampling validation remain open. |
| `reference_fcp` | Reference-only | Reference-algorithm reproduction and historical real-engine parity, with the later energy-handoff correction checked separately. Not accepted for relax/NEB/MD. |

The canonical configuration boundary is in the
[configuration reference](../user/CONFIG_REFERENCE.md#211-constant-potential-evaluation).
The current feature matrix remains the support-status entry point.

## Historical relaxation and NEB acceptance recovered into ATST

The [compact evidence records](data/ATST_CP_WORKFLOW_20260921/README.md)
were retained in the parent ABACUS project's acceptance directory but were
missing from ATST's current status summary. They establish more than a mock
workflow or a single successful SCF:

- Relaxation restart job **1425120** completed an interruption and restart,
  advanced two nuclear steps, and reached fmax **0.025953 eV/Angstrom** against
  the original **0.05** criterion. The final potential residual was
  **1.03e-8 eV** or less. The first resumed electron count equaled the saved
  evaluated count.
- Serial NEB job **1425552** produced **176 frames / 44 complete bands**.
  Final fmax was **0.070402 eV/Angstrom**, below the original **0.08** criterion;
  the maximum final potential residual was **6.75e-9 eV**. The fixed Pt
  substrate stayed fixed to floating-point precision. A read-only trajectory
  check on 2026-10-20 reconfirmed its identity, frame count, fixed-atom geometry
  and final potential residual.
- NEB restart job **1426458** resumed each interior image from its own
  evaluated electron count and advanced six FIRE steps. Final fmax was
  **0.077939 eV/Angstrom** and the maximum potential residual was
  **1.14e-5 eV**, meeting the original **0.08 / 1e-4** criteria.

These runs used the historical implementation integrated in `d30747b`, which
is an ancestor of 2.3.2. They were not rerun here with the current package.
The independent-backend migration and energy-handoff repairs occurred later;
their [release evidence](../releases/RELEASE_NOTES_2.3.1.md) and current
[2.3.2 installed-consumer verification](../releases/RELEASE_NOTES_2.3.2.md)
are separate evidence. The final-source MD GPU run exercises the current
compensated solver and moving geometry, but does not replace a current-version
full relax/NEB rerun.

The historical optimizer restarts rebuild BFGS/FIRE rather than restoring
optimizer state. Both interrupted and uninterrupted runs met the original
criteria, with different paths/final coordinates. The small two-interior-image
NEB fixture does not establish a converged meV-scale barrier.

## Numerical and application boundaries

The [Pt example](../../examples/19_constant_potential_Pt/README.md) records
compensated chemical-potential derivative and reconverged-force checks, plus
PW numerical and charge-confinement limitations. The
[runtime joint report](ATST_CP_RUNTIME_JOINT_VALIDATION_2026-09-21.md),
[reference parity report](ATST_CP_PARITY_2026-10-20.md) and
[MD validation report](ATST_CP_MD_VALIDATION_2026-10-20.md) retain their dated
scope; this reassessment does not retroactively change their conclusions.

Supported calculations still need an appropriate Hamiltonian, charge
confinement, k-mesh/basis/cutoff convergence and observable-specific sampling.
The custom compensated chemical potential has no automatic SHE/RHE mapping.
The short real MD runs did not resolve clean second-order timestep convergence;
longer real-backend validation should examine integration error, electronic
convergence and confinement together. Stochastic restarts do not restore RNG
state, and output-tail inconsistencies are rejected rather than repaired.

The remaining targeted verification work is to refresh the historical
relax/serial-NEB case with the current installed backend, exercise CP-specific
image-parallel execution, and extend real-backend MD integration checks. This
list does not imply support for NPT, variable-cell CP, implicit solvent or other
configurations rejected by the current schema/factory.

## Downstream and compatibility boundaries

Historical SAI acceptance included a packaged prepare/run/summary chain for
relaxation and serial NEB. That establishes the recorded tool handoff, not a
current deployed Paimon/Toolbox release or platform Flow acceptance. Downstream
release status belongs to its owner and does not determine ATST's software
support label.

The persisted `profile_status` strings remain unchanged for facts-v1 reader
compatibility. They are historical labels, not a live software-maturity
classification. Current user documentation and the feature matrix carry that
classification; the actual SCF/CP convergence and identity fields continue to
determine whether a saved evaluation is usable.
