# Constant-potential MD implementation and validation — 2026-10-20

## Delivery scope

The development branch adds fixed-cell ASE MD to the existing compensated-gate
constant-potential calculator. It reuses ASE VelocityVerlet, Bussi, Langevin and
NVTBerendsen. Native ABACUS MD, NPT/variable cell, reference_fcp and target lists
are rejected. Each force comes from a converged SCF and electron-number solve.
The previous complete electron count is an initial guess only.

This is an unreleased extension based on main `d41eba7`. The installed dependency
remains `abacuslite>=1.2.1,<2`; no backend source change or ABACUS PR #8136 is
required. Published ATST 2.3.1 does not contain these changes.

The lasting user contract is in the [configuration reference](../user/CONFIG_REFERENCE.md#experimental-constant-potential-md).
The [Pt materializer](../../examples/19_constant_potential_Pt/materialize_md.py)
requires an explicit target for its own Hamiltonian, fixes eight Pt atoms and
leaves four mobile. It uses LCAO/cusolver, a 2 Ry block and explicit gate/dipole
settings. It is a bounded workflow fixture, not a converged electrode model.

## Offline verification and independent review

- Full offline suite: **1450 passed, 37 skipped, 1 warning**. Skips are existing
  opt-in MPI, external DeePMD, built-distribution and unavailable-fixture checks.
- Analytic electronic model: geometry-dependent stationary electron count,
  minimized grand-potential forces and finite differences; real production CP
  iterations and ASE integration; constraints applied once; no extra solves
  from output reads; all four integrators; timestep refinement and deterministic
  split/restart comparisons.
- Persistence: step-zero and every-complete-step atomic checkpoints; failed SCF
  and CP solves retain the previous complete state; finite zero momenta survive
  restart; interval-five output resumes on global cadence after step seven.
- Offline postprocessing validates saved CP facts and uses saved step/time,
  without attaching an electronic calculator.
- Independent review found and repaired two issues: a direct native-runner
  bypass and incomplete cross-validation of restart output tails. The final
  review accepted unconditional native rejection and nuclear-state fingerprints
  plus saved-energy checks. Restart preflight was independently checked with
  calculator construction and CP evaluation forbidden.

The source review snapshots and raw test logs are retained under
`/home/james/scratch/atst-abacuslite-release-20261020/`:
`cp-md-review-fixed/manifest.json`, `cp-md-final.log`, `cp-md-final.xml`, and
`cp-md-runtime-review-green1.log`. Source hashes bind the dirty workstate rather
than treating the baseline HEAD as the changed implementation.

## GPU verification

Job **1872958** completed with exit `0:0` in 19 min 24 s on SAI `4V100`,
node `4v100n14`, one V100, eight OpenMP threads. It used ABACUS
`LTSv3.10.1-sm70-auto`, Python 3.10.19, external abacuslite 1.2.1, ASE 3.28.0
and NumPy 2.2.6.
The source/runner bundle and complete evidence archive are retained under
`/org/pku-jianghong/liuzhaoqing/abacus-cp-development/20261020-cp-md/`;
local evidence is under the corresponding `cp-md-gpu/` scratch directory.

The target was measured from the exact initial fixed-N calculation at N=217:
`mu = 15.04601576397786 eV`, with N0=216. Both gate and dipole derivatives
were read from that calculation's charge density. There is no SHE/RHE mapping.

| Observation | Result |
| --- | --- |
| Gate-facing Pt x-force at the initial geometry | 1.6568329131 eV/Å |
| Centered grand-potential difference, h=0.01 Å | 1.6574422472 eV/Å |
| Force minus finite difference | −0.0006093341 eV/Å |
| Langevin NVT, 2 × 0.25 fs, target 100 K | Completed; maximum CP residual 2.14×10⁻⁶ eV |
| VelocityVerlet, 2 × 0.25 fs | Energy range 5.0753×10⁻⁶ eV |
| VelocityVerlet, 4 × 0.125 fs | Energy range 5.0694×10⁻⁶ eV |
| Split 1+1 versus continuous 2-step Verlet | Maximum position difference 3.11×10⁻¹¹ Å; momentum difference 1.26×10⁻⁷ in ASE units |
| Looser CP tolerance, 10⁻⁴ versus 10⁻⁵ eV, same 0.25 fs | Energy range 5.0653×10⁻⁶ eV; actual residual still about 2.41×10⁻⁶ eV |
| Charge in the gate/block interval, x/Lx∈[0.65,0.75), across 38 SCFs | 0.0001535–0.0001636 electrons |

The displacement tests reconverged electrons at each geometry. Halving the MD
timestep did not produce a clean second-order improvement in this very short
real-backend run: the roughly 5 micro-eV variation cannot isolate integration
error from backend force/energy precision. Likewise, loosening the electronic
stopping threshold did not produce a comparably looser actual residual. These
are observations, not claims of long-time convergence or tolerance independence.
The NVT energy range is not an energy-conservation criterion.

This first job used the pre-final-review execution snapshot. The electronic
solver and integrators are unchanged, but final restart identity validation was
subsequently tightened. Job **1873489** reran the final reviewed source on
`4v100n17` and completed with exit `0:0` in 7 min 4 s:

- One 0.25 fs Langevin step and two 0.25 fs VelocityVerlet steps completed.
- A separate 1+1-step Verlet continuation passed the final nuclear-state and
  saved-energy checks. Its final position differed from the continuous run by
  at most **7.41×10⁻¹² Å**, with momentum difference **1.87×10⁻⁸** in ASE units.
- Maximum CP residual across these runs was **2.46×10⁻⁶ eV**, below 10⁻⁵ eV.
- Saved nuclear frames keep x between 4.99997 and 9.53087 Å, outside the
  gate/block interval beginning near 13 Å.
- Every completed diagnostic row contains the final nuclear-state fingerprint.
  The complete source manifest was compared with the delivered source: no drift.

Final-job evidence is retained under the sibling site directory
`20261020-cp-md-final/` (`evidence-1873489.tar.gz`, scheduler logs and input
bundle); local summary, source manifest, script and logs are in `cp-md-gpu-final/`.
The exact uploaded source bundle SHA-256 is
`780491be66ae0701ab7b858ee40186bfe30ef2e10e4b4be26f05104e1792dd0b`.
Both jobs used ABACUS executable SHA-256
`c29ff3096b9edf455e76e36e7055ef3937c56194bc0235ae9df5e8cfb8c3c303`.
The full charge-density archives remain on SAI; local JSON/logs retain the
reported observations without requiring a copy of all intermediate SCFs.

## Scientific and restart limits

The diagnostic energy is the force-consistent grand potential plus nuclear
constraint energy, with nuclear kinetic energy added for the unthermostatted
energy check. Raw DFT energy plus kinetic energy is not conserved under electron
exchange. Berendsen is not exact canonical sampling. Stochastic thermostat
state/RNG is not checkpointed, so restart is physical state continuation rather
than guaranteed bitwise stochastic continuation.

The checkpoint is replaced atomically per complete MD step. Trajectory, JSONL
and log files are separate writes: an I/O interruption between writes may leave
inconsistent tails. Restart detects and rejects that inconsistency; it does not
silently truncate or repair output. Electronic failure inside an integration
step preserves the preceding complete checkpoint.

No merge, package-version bump, publication or upstream PR comment is part of
this implementation delivery. A later release must add its own release notes
and verify the exact integrated revision.
