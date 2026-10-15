"""Fixed-cell input validation and ASE NEB runtime preparation."""

from __future__ import annotations

from copy import deepcopy
import warnings

import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator

from atst_tools.calculators.constant_potential import clear_constant_potential_facts
from atst_tools.utils.neb_endpoints import ENDPOINT_RESULT_KEY


CELL_MATCH_ATOL = 1e-5  # Absolute Angstrom tolerance; no relative term.
CELL_PREPARATION_KEY = "atst_cell_preparation"
# ASE 3.28 checks periodic vector components at this bound. This is an
# adaptation trigger, not a replacement for ASE's native constructor check.
_ASE_NEB_CELL_ATOL = 1e-8


def copy_with_cached_results(atoms):
    """Copy an image and its valid frozen results without evaluating a backend.

    Args:
        atoms: Source ASE image, possibly carrying a SinglePointCalculator.

    Returns:
        An independent image with deeply copied metadata and frozen results.
        Live calculators are not shared with the source image.
    """
    image = atoms.copy()
    image.info = deepcopy(atoms.info)
    if isinstance(atoms.calc, SinglePointCalculator) and not atoms.calc.check_state(atoms):
        image.calc = SinglePointCalculator(image, **deepcopy(atoms.calc.results))
    return image


def prepare_fixed_cell_images(images, *, context="Fixed-cell input", runtime=False, world=None, reference=None):
    """Validate boundaries and prepare private image copies against image zero.

    Args:
        images: Nonempty sequence whose first image defines the reference.
        context: English diagnostic label identifying the workflow.
        runtime: Also adapt periodic vectors that exceed ASE's NEB bound.
            Input preparation otherwise preserves all cell values within
            ``CELL_MATCH_ATOL``. Cartesian positions are never rescaled.
        world: Existing communicator used only to restrict diagnostics to root.
            This helper never creates or initializes a communicator.
        reference: Optional common input reference for segmented preprocessing.
            Defaults to the first image and never changes within this call.

    Returns:
        Private images. Changed images carry their original boundaries and
        the actual runtime boundaries in ``CELL_PREPARATION_KEY``; cached
        results and constant-potential identity are discarded on any change.

    Raises:
        ValueError: If the sequence is empty, coordinates/cells are non-finite,
            or reference periodic vectors are missing or linearly dependent.
    """
    if not images:
        raise ValueError(f"{context} requires at least one image")
    reference = images[0] if reference is None else reference
    reference_cell = reference.cell.array.copy()
    reference_pbc = reference.pbc.copy()
    if not np.isfinite(reference_cell).all() or not np.isfinite(reference.positions).all():
        raise ValueError(f"{context} reference image 0 has non-finite cell or positions")
    for index, atoms in enumerate(images):
        if not np.isfinite(atoms.cell.array).all() or not np.isfinite(atoms.positions).all():
            raise ValueError(f"{context} image {index} has non-finite cell or positions")
    periodic_vectors = reference_cell[reference_pbc]
    if len(periodic_vectors) and np.linalg.matrix_rank(periodic_vectors) < len(periodic_vectors):
        raise ValueError(f"{context} reference image 0 has invalid periodic cell vectors")
    # AutoNEB may make any interior image the first image of a later sub-band.
    # Bound the span after input overrides, rather than only differences to
    # image zero: two individually acceptable offsets can exceed ASE's guard.
    prospective_cells = [
        reference_cell if (
            np.max(np.abs(atoms.cell.array - reference_cell)) > CELL_MATCH_ATOL
            or not np.array_equal(atoms.pbc, reference_pbc)
        ) else atoms.cell.array
        for atoms in images
    ]
    runtime_requires_common_cell = runtime and np.any(
        np.ptp(np.asarray([reference_cell, *prospective_cells])[:, reference_pbc], axis=0)
        > _ASE_NEB_CELL_ATOL
    )
    prepared = []
    for index, atoms in enumerate(images):
        image = copy_with_cached_results(atoms)
        original_cell = image.cell.array.copy()
        original_pbc = image.pbc.copy()
        max_delta = float(np.max(np.abs(original_cell - reference_cell)))
        pbc_mismatch = not np.array_equal(original_pbc, reference_pbc)
        cell_mismatch = max_delta > CELL_MATCH_ATOL
        reasons = []
        if cell_mismatch or pbc_mismatch:
            if world is None or world.rank == 0:
                warnings.warn(
                    f"{context} image {index}: adopting reference image 0 cell/PBC; "
                    f"max cell delta={max_delta:.9g} Angstrom, atol={CELL_MATCH_ATOL:g}, "
                    f"PBC={original_pbc.tolist()} -> {reference_pbc.tolist()}; "
                    "Cartesian positions are unchanged.",
                    UserWarning, stacklevel=2,
                )
            image.set_cell(reference_cell, scale_atoms=False)
            image.set_pbc(reference_pbc)
            reasons.append("input_boundary_override")
        elif runtime_requires_common_cell and not np.array_equal(original_cell[reference_pbc], periodic_vectors):
            cell = original_cell.copy()
            cell[reference_pbc] = periodic_vectors
            image.set_cell(cell, scale_atoms=False)
            reasons.append("ase_neb_fixed_cell")
        if reasons:
            image.calc = None
            image.info.pop(ENDPOINT_RESULT_KEY, None)
            clear_constant_potential_facts(image)
            record = {
                "reference_image": 0,
                "reference_cell": reference_cell.tolist(),
                "reference_pbc": reference_pbc.tolist(),
                "original_cell": original_cell.tolist(),
                "original_pbc": original_pbc.tolist(),
                "actual_cell": image.cell.array.tolist(),
                "actual_pbc": image.pbc.tolist(),
                "max_cell_delta_angstrom": max_delta,
                "atol_angstrom": CELL_MATCH_ATOL,
                "reasons": reasons,
            }
            image.info[CELL_PREPARATION_KEY] = [
                *image.info.get(CELL_PREPARATION_KEY, []), record,
            ]
        prepared.append(image)
    return prepared


def fixed_cell_metadata(images):
    """Return durable reference and adaptation facts for a prepared band.

    Args:
        images: Prepared images in workflow order.

    Returns:
        JSON-compatible reference boundaries and per-image adjustments.
    """
    return {
        "reference_image": 0,
        "reference_cell": images[0].cell.array.tolist(),
        "reference_pbc": images[0].pbc.tolist(),
        "atol_angstrom": CELL_MATCH_ATOL,
        "adjustments": [
            {"image": index, **record}
            for index, image in enumerate(images)
            for record in image.info.get(CELL_PREPARATION_KEY, [])
        ],
    }
