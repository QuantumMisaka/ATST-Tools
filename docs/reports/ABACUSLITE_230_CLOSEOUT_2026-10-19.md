# ATST 2.3.0 integration closeout

Date: 2026-10-19. Scope: the external abacuslite integration, installed
distribution acceptance, and disposition of local development branches.
ATST 2.3.0 is published on PyPI (2026-10-19); publication evidence is recorded in the
[release notes](../releases/RELEASE_NOTES_2.3.0.md).

## Update 2026-10-19 (user decision)

The release candidate now targets the published `abacuslite==1.2.0`
(wheel sha256 `3e5e245c0a03e8182c2c726d7ccad1c800c20250641abd4b81a8e22e892bc4a0`,
PyPI upload 2026-10-19T03:41:48Z). The supported dependency range is
`abacuslite>=1.2.0,<2`, and 2.3.0 installed-distribution plus installed-consumer
acceptance passed against 1.2.0 before publication. The `1.1.0` acceptance
below stays as the historical record for the earlier candidate; it is not
relabeled and its artifacts are not reused for the 1.2.0 candidate.

## Historical 1.1.0 candidate baseline

The historical runtime candidate was ATST `49c8626a9d85858f7022fc7a732cbb01d871298a`
with published `abacuslite==1.1.0` (tag source
`94a26e4c7f31ed817eac36c871ff74a7b7abe418`). That candidate used the dependency range
`abacuslite>=1.1.0,<2`; the published release requires 1.2.0 as stated above. ATST owns workflow, launcher and result
contracts; the independent package owns calculator and parser behavior.
Historical imports are forwarding modules, with no vendored runtime fallback.

Historical ruling (superseded by the 1.2.0 release decision above): retain the
published 1.1.0 floor for that candidate. The plugin's
post-tag PP/ORB directory defaults, explicit atomic-mass preservation and
configured restart-input repairs are not included in the 1.1.0 wheel. A future
ATST requirement for those behaviors must follow a plugin release and targeted
consumer acceptance, rather than copying the fixes into ATST or relying on an
unpublished plugin checkout.

## Historical 1.1.0 candidate acceptance

- Main-branch [GitHub run 39145038570](https://github.com/QuantumMisaka/atst-tools/actions/runs/39145038570)
  passes four Python/ASE lanes, documentation governance and the installed
  external consumer lane at `49c8626`.
- Local source unit tests pass with five environment-gated skips when
  `PYTHONPATH=src` selects this checkout. The shared `atst-dev` environment's
  default editable import points to another checkout; its default import is
  not evidence for this worktree.
- A wheel built from the source distribution at `49c8626` is retained at
  `/home/james/scratch/atst-230-closeout-20261019.NjC3Jb/dist/`.
  Wheel SHA256: `a127b3972458abfdda33503e16c2780b0ecdb98e53136d22fe4dfa97fdbc7301`.
  Sdist SHA256: `277c1ee5ca0b44dd3317d10963a61030341fa5550cc1519c4eb7dd76ec67f793`.
- A fresh Python 3.10 virtual environment with ASE 3.28.0 and the published
  plugin passes `pip check`, strict wheel/sdist metadata checks and the full
  installed-consumer suite: 1419 collected, 1386 passed, 33 environment-gated
  skips, zero failures or deselections. Both distribution artifacts were
  supplied to the suite; the consumer report verifier passes. Reports and
  raw test logs are retained under
  `/home/james/scratch/atst-230-closeout-20261019.NjC3Jb/`.
  Python 3.10's bundled pip 23.0.1 omits a bare local wheel archive hash;
  installation used a `file:` wheel URL with an explicit SHA256 fragment,
  which satisfies the existing provenance gate without changing it.
- `scripts/verify_wheel_api.py --wheel ...` passes its separate clean-install
  CLI/API dry runs, stored-data TS validation and public API example checks.
- Historical 1.1.0-candidate GPU acceptance passed as SAI job 1848462 on `4v100n09`
  (4V100, two GPUs): `COMPLETED`, exit `0:0`, elapsed 11m49s. All 18 stages
  and 13 real-engine cases pass; the retained-archive audit returns
  `pass: true`, `errors: []`. The develop engine source is `4b92455` and
  the legacy-reader lane uses ABACUS LTS 3.10.1. Evidence is retained in
  `runtime/` beneath the acceptance directory above. Result archive SHA256:
  `4565d61aa42e7ed1878b44d754728bc27e7d0f5148b129599fb3392a5e7bb3e7`.
  The earlier SAI job 1835021 remains historical candidate evidence and is
  not relabeled as this 2.3.0 wheel's result.

## Branch and worktree disposition

| Branch | Evidence | Disposition |
| --- | --- | --- |
| `integrate/abacuslite-atst229-20261017` at `07493db` | Nine commits differ by history from main. Native MD, CP SCF markers and gate-cube handling are already identical to main. The old MD cell/PBC reader repair exists in published abacuslite 1.1.0. Remaining adapter differences are the superseded vendored fallback, provenance and older command probing. Independent review found no implementation to transplant. | Superseded historical candidate; do not merge it wholesale. Retain the branch and Codex-managed worktree at `/home/james/.codex/worktrees/5738/atst-tools`, including ignored scientific artifacts. |
| `fix/ase-ts-dual-interface-20261013` at `3163c06` | Commit is an ancestor of main; tracked and untracked state is clean, but ignored scientific artifacts remain. | Completed historical worktree, already archived under `/home/james/scratch/atst-tools-closeout-20261015-153534/atst-tools-worktree`; retain its branch and artifacts. |
| `docs/gpu-p6-inner-mpi-20261016`, `fix/review-nine-findings-20261016`, `release/2.2.9-p6-closeout-20261017` | All are ancestors of main; none has a separate checked-out worktree. | Renamed locally to `archive/gpu-p6-inner-mpi-20261016`, `archive/review-nine-findings-20261016` and `archive/2.2.9-p6-closeout-20261017`. Every original commit is retained; the reversible mapping is in `branch-archive.json` beside the acceptance artifacts. |

At the original closeout, the runtime baseline was synchronized to `49c8626`; it added only
documentation. The parent
`app-tools` index then recorded `45d8c78`; the parent owner must integrate the
chosen submodule revision with its own delivery changes. No parent gitlink
commit is implied by this report.

## Publication workspace finalization (2026-10-19)

PRs #39, #40 and #41 are merged. Their branch file trees were compared with
merge commits `fa1ba07`, `b8c3218` and `2465e24`, respectively, with no differences.
Local recovery branches now preserve the three heads:

- `archive/230-integration-closeout-20261019` retains `d7f8971`.
- `archive/230-abacuslite-120-release-20261019` retains `e58900f`.
- `archive/230-publication-record-20261019` retains `f287091`.

The two historical worktrees listed above remain in place with their scientific
artifacts. At that checkpoint, the parent-managed submodule checkout remained at `45d8c78`;
publication finalization uses an isolated worktree and does not update the
parent gitlink. Local `main` has been fast-forwarded to the published-source
history. The reversible branch mapping and verification artifacts are retained
under `/home/james/scratch/atst-230-finalize-20261019-za20o8qf/`.

The published-wheel consumer recheck with abacuslite 1.2.0 passed all 1386
executed tests (33 environment-gated skips, no failures or deselections), with
both PyPI distribution artifacts supplied. Package imports resolve to the
isolated environment, `BACKEND_SOURCE=external`, and `pip check` passes.
This recheck does not add GPU or scientific-convergence evidence.

## Integration and remaining delivery work

| Work | Owner and completion boundary |
| --- | --- |
| Plugin consumer CI reference | [abacuslite PR #5](https://github.com/QuantumMisaka/abacuslite/pull/5) replaces the deleted ATST integration branch with main and pins PR checks to `49c8626`. The follow-up documentation correction distinguishes proposed ABACUS adoption from merged ATST adoption. All 16 checks passed at `6409df6`; merged as `87f0960` on 2026-10-19 (Asia/Shanghai). |
| ATST 2.3.0 publication | Completed 2026-10-19: tag `v2.3.0` at `b8c3218`, publish run [39196721131](https://github.com/QuantumMisaka/atst-tools/actions/runs/39196721131) success, PyPI wheel `80a5e242…5729` and sdist `af317351…3bf2`; fresh-venv install reports 2.3.0 with `BACKEND_SOURCE=external`. See the [release notes](../releases/RELEASE_NOTES_2.3.0.md). |
| Upstream engine adoption | [ABACUS PR #8175](https://github.com/deepmodeling/abacus-develop/pull/8175) is the external-only engine/harness adoption. After upstream merge, the plugin owner updates its reviewed harness reference and reruns its engine gate. |
| GPU P6 / Paimon delivery | Parent ABACUS/MLIP-Agent owners retain the [active platform plan](../superpowers/plans/2026-09-20-atst-gpu-node-tuning-plan.md). Package acceptance does not close SIF, platform invocation or production publication. |

The runtime harness covers execution and artifact handoff for Relax, ASE MD
initial/restart, native MD, serial/parallel NEB, CP baseline/point, vibration,
legacy-reader Relax, Sella, CCQN and descent IRC. Dimer and scientific
convergence studies are excluded. One-step optimizer completion does not
establish a converged structure, transition state or reaction path.

## Downstream status verified on 2026-10-20

The parent index now pins ATST `e427a70` and abacuslite `e8acc4c` (2.3.0 /
1.2.0 source lines). Parent P6 code is merged. The parent's dated
`toolbox/ABACUS/docs/reports/2026-10-19-abacus-1.2.76-final-platform-acceptance.md`
records Toolbox 1.2.76 deployed to dev/prod with Agent and Flow GaN
relax/band/DOS acceptance. Its SIF contains abacuslite 1.2.0 and ATST is
delivered with the Toolbox. These facts supersede the earlier parent-pin and
general deployment blockers; they do not close P6's 1/4-card NEB comparison.
The 1.2.1 / 2.3.1 package closeout does not update those downstream pins or SIF.
