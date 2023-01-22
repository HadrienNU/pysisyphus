# [1] https://doi.org/10.1063/1.3463717
#     Comparisons of classical and Wigner sampling of transition state
#     energy levels for quasiclassical trajectory chemical dynamics simulations
#     Sun, Hase, 2010
# [2] https://doi.org/10.1002/qua.25049
#     Effects of different initial condition samplings on photodynamics and
#     spectrum of pyrrole
#     Barbatti, Sen, 2015

import functools
from typing import Callable

import numpy as np
from numpy.typing import NDArray

from pysisyphus.constants import AMU2AU, AMU2KG, BOHR2M, C, HBAR, P_AU
from pysisyphus.helpers_pure import eigval_to_wavenumber
from pysisyphus.Geometry import Geometry, get_trans_rot_projector


def normal_mode_reduced_masses(masses_rep, normal_modes):
    return 1 / (np.square(normal_modes) / masses_rep[:, None]).sum(axis=0)


@functools.singledispatch
def get_wigner_sampler(
    coords3d: NDArray,
    masses: NDArray,
    hessian: NDArray,
    nu_thresh: float = 5.0,
    stddevs: float = 6.0,
) -> Callable:
    assert coords3d.shape == (len(masses), 3)
    assert hessian.shape == (coords3d.size, coords3d.size)

    # Projector to remove translation & rotation
    Proj = get_trans_rot_projector(coords3d, masses, full=True)
    masses_rep = np.repeat(masses, 3)
    mm_sqrt = np.sqrt(masses_rep)
    M = np.diag(mm_sqrt)
    M_inv = np.diag(1 / mm_sqrt)
    PM = Proj @ M_inv

    # Diagonalize projected, mass-weighted Hessian.
    w, v = np.linalg.eigh(PM @ hessian @ PM.T)
    nus = eigval_to_wavenumber(w)
    small_nu_mask = np.abs(nus) < nu_thresh
    w = w[~small_nu_mask]
    v = v[:, ~small_nu_mask]
    nus = nus[~small_nu_mask]
    assert (nus >= nu_thresh).all(), "Imaginary wavenumbers are not yet handled!"
    n = len(nus)  # Number of non-zero wavenumbers

    # Convert wavenumbers (cm⁻¹) to angular frequency in 1/s
    ang_freqs = 2 * np.pi * nus * 100 * C

    # Reduced masses in amu and then in kg
    mus_amu = normal_mode_reduced_masses(masses_rep, v)
    mus = mus_amu * AMU2KG

    # Standard deviations. Eq. (3) in [2].
    sigma2_q = HBAR / (2 * mus * ang_freqs)
    sigma2_p = (HBAR * mus * ang_freqs) / 2
    # Sigmas are now in atomic units!
    sigma2_q = sigma2_q / (BOHR2M**2)
    sigma2_p = sigma2_p / (P_AU**2)
    q_quot = 1 / (2 * sigma2_q)
    p_quot = 1 / (2 * sigma2_p)

    span = 2 * stddevs

    def sampler():
        Qs = np.zeros(n)
        Ps = np.zeros(n)
        for i in range(n):
            qi_quot = q_quot[i]
            pi_quot = p_quot[i]
            while True:
                q, p, ref = np.random.random_sample(3)
                # Map q and p from [0, 1) onto chosen interval [-stddevs, stddevs)
                q = q * span - stddevs
                p = p * span - stddevs
                # Wigner function. Eq. (2) in [2].
                p_wig = (
                    1
                    / np.pi
                    * np.exp(-(q**2) * qi_quot)
                    * np.exp(-(p**2) * pi_quot)
                )
                if p_wig >= ref:
                    break
            Qs[i] = q
            Ps[i] = p

        # Convert to Cartesian coordinates
        displ = M_inv @ v @ Qs
        velocities = M @ v @ Ps

        # The COM remains unaffected, as we displace along vibrations,
        # not translations.
        displ_coords3d = coords3d.copy() + displ.reshape(-1, 3)

        # dp is still in atomic units with electron_mass as mass unit
        velocities /= masses_rep * AMU2AU
        P = get_trans_rot_projector(displ_coords3d, masses)
        velocities = P.dot(velocities)
        return displ_coords3d, velocities

    return sampler


@get_wigner_sampler.register
def _(geom: Geometry):
    return get_wigner_sampler(geom.coords3d, geom.masses, geom.hessian)
