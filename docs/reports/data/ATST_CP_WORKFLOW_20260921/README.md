# Historical constant-potential optimization evidence

These compact records preserve the September 2026 relax/serial-NEB acceptance
that was recorded in the parent ABACUS project. They support the
[2026-10-20 support assessment](../../ATST_CP_SUPPORT_ASSESSMENT_2026-10-20.md).
They are historical observations, not new runs of ATST 2.3.2.

- [neb-summary.json](neb-summary.json) projects the completed job 1425552 observation. On
  2026-10-20, the retained `neb.traj` was read to recheck its SHA256, 176 frames,
  44 bands, fixed-Pt displacement and final CP residual. The displacement here
  is the maximum Euclidean norm (4.97e-16 Angstrom); the original observation
  reported the maximum component (4.44e-16 Angstrom). NEB fmax and grand
  energies are retained from the historical observation, not recomputed by a
  new optimizer or electronic calculation.
- [restart-relax-controller-review.json](restart-relax-controller-review.json) is an unchanged copy of the original
  controller check. Job 1425101 rebuilt an already converged state; job 1425120
  also continued two nuclear steps. Only the latter establishes advancing
  continuation.
- [restart-neb-controller-review.json](restart-neb-controller-review.json) is an unchanged copy of the original
  controller check of job 1426458. Each interior image reused its own evaluated
  electron count; resumed optimization met the original convergence criteria.

Original records and trajectories remain under
`/home/james/scratch/constant-potential-sdd-20260920/paimon-sai/`.
Full electronic outputs remain in the historical SAI job directories. The
parent development record is
`toolbox/ABACUS/docs/superpowers/plans/2026-09-20-abacus-constant-potential-development-plan.md`,
entries for jobs 1425120, 1425552 and 1426458. The ATST implementation was
integrated in `d30747b6c7f826d7f53e6fc31e9534ad9e2c36ed`, an ancestor of the
2.3.2 release. Later source/backend changes and current verification are
accounted for separately in the assessment.

SHA256 of the unchanged restart records:

| File | SHA256 |
| --- | --- |
| `restart-relax-controller-review.json` | `fb91ee05428b3625ed88ee53425dba497db71b67dce08409730620f0296a4bac` |
| `restart-neb-controller-review.json` | `bbd2ac8a6abd38dfafc5cd6a2980c1031b00bdb94c1416435374b33d68f139a0` |
