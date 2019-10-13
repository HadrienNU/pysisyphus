#!/usr/bin/env python3

# [1] https://doi.org/10.1016/S0022-2860(84)87198-7
#     Pulay, 1984
# [2] https://pubs.rsc.org/en/content/articlehtml/2002/cp/b108658h
#     Farkas, Schlegel, 2002

from collections import namedtuple
import logging

import autograd.numpy as anp
from autograd import grad
import numpy as np
from scipy.optimize import minimize


COS_CUTOFFS = {
    2: 0.97,
    3: 0.84,
    2: 0.80,
    3: 0.75,
    4: 0.71,
    5: 0.67,
    6: 0.62,
    7: 0.56,
    8: 0.49,
    9: 0.41,
}
DIISResult = namedtuple("DIISResult",
                         "coeffs coords forces N"
)
logger = logging.getLogger("optimizer")


def log(msg):
    logger.debug(msg)


def valid_diis_direction(diis_step, ref_step, use):
    ref_direction = ref_step / np.linalg.norm(ref_step)
    diis_direction = diis_step / np.linalg.norm(diis_step)
    cos = diis_direction @ ref_direction
    return (cos >= COS_CUTOFFS[use]) and (cos >= 0)


def from_coeffs(vec, coeffs):
    return np.sum(coeffs[:,None] * vec[::-1][:len(coeffs)], axis=0)


def diis_result(coeffs, coords, forces):
    diis_coords = from_coeffs(coords, coeffs)
    diis_forces = from_coeffs(forces, coeffs)
    diis_result = DIISResult(
                        coeffs=coeffs,
                        coords=diis_coords,
                        forces=diis_forces,
                        N=len(coeffs),
    )
    log(f"\tUsed {len(coeffs)} error vectors for DIIS.")
    log("")
    return diis_result


def gdiis(err_vecs, coords, forces, ref_step, max_vecs=5):

    # Scale error vectors so the smallest norm is 1
    norms = np.linalg.norm(err_vecs, axis=1)
    err_vecs = err_vecs / norms.min()

    valid_coeffs = None
    for use in range(2, min(max_vecs, len(err_vecs))+1):
        log(f"Trying GDIIS with {use} previous cycles.")
        use_vecs = np.array(err_vecs[::-1][:use])

        A = np.einsum("ij,kj->ik", use_vecs, use_vecs)
        coeffs = np.linalg.solve(A, np.ones(use))
        # Scale coeffs so that their sum equals 1
        coeffs_norm = np.linalg.norm(coeffs)
        valid_coeffs_norm = coeffs_norm <= 1e8
        log(f"\tError vectors are linearly independent: {valid_coeffs_norm}")
        coeffs /= np.sum(coeffs)
        coeffs_str = np.array2string(coeffs, precision=4)
        log(f"\tGDIIS coefficients: {coeffs_str}")

        # Uncomment these lines and break here to only do the basic check
        # for linear dependency above.
        # valid_coeffs = coeffs
        # break

        # Check degree of extra- and interpolation.
        pos_sum = abs(coeffs[coeffs > 0].sum())
        neg_sum = abs(coeffs[coeffs < 0].sum())
        valid_sums = (pos_sum <= 15) and (neg_sum <= 15)
        log(f"\tSum of positive coefficients: {pos_sum:.2f}")
        log(f"\tSum of negative coefficients: {neg_sum:.2f}")
        log(f"\tSums are valid: {valid_sums}")

        # Calculate GDIIS step for comparison to the reference step
        diis_coords = from_coeffs(coords, coeffs)
        diis_step = diis_coords - coords[-1]
        valid_length = np.linalg.norm(diis_step) <= (10 * np.linalg.norm(ref_step))
        log(f"\tGDIIS step has valid length: {valid_length}")

        # Compare directions of GDIIS- and reference step
        valid_direction = valid_diis_direction(diis_step, ref_step, use)
        log(f"\tGDIIS step has valid direction: {valid_direction}")

        gdiis_valid = (valid_sums
                       and valid_coeffs_norm
                       and valid_direction
                       and valid_length
        )
        log(f"\tGDIIS step is valid: {gdiis_valid}")
        if not gdiis_valid:
            break
        # Update valid DIIS coefficients
        valid_coeffs = coeffs
        log("")

    if valid_coeffs is None:
        return None

    return diis_result(valid_coeffs, coords, forces)


def gediis(coords, energies, forces, max_vecs=10):
    use = min(len(coords), max_vecs)

    R = coords[::-1][:use]
    E = energies[::-1][:use]
    f = forces[::-1][:use]
    # Precompute values so they can be reused in fun()
    Rifi = np.einsum("ik,ik->i", R, f)
    Rjfi = np.einsum("jk,ik->ji", R, f)

    def x2c(x):
        return x**2 / (x**2).sum()

    # def fun(xs):
        # """Naive implementation with loops."""
        # cs = x2c(xs)
        # first = (cs*E).sum()
        # sec = 0.
        # for i, ci in enumerate(cs):
            # for j, cj in enumerate(cs):
                # sec += ci * cj * (f[j] - f[i]) @ (R[i] - R[j])
        # return first - 1/2 * sec

    # def fun(xs):
        # """Recalculation of all values in every call."""
        # cs = x2c(xs)
        # return anp.sum(cs*E) - anp.einsum("i,j,jk,ik", cs, cs, R, f) + anp.einsum("i,ij,ij", cs, R, f)

    def fun(xs):
        """Usage of precomputed values."""
        cs = x2c(xs)
        return anp.sum(cs*E) - anp.sum(anp.outer(cs, cs)*Rjfi) + (cs * Rifi).sum()

    jac = grad(fun)

    x0 = np.ones(use)
    res = minimize(fun, x0=x0, jac=jac)
    # print(res)

    coeffs = None
    if res.success:
        coeffs = x2c(res.x)
    return diis_result(coeffs, coords, forces)
