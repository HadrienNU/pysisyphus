#!/usr/bin/env python3

from collections import namedtuple
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.TablePrinter import TablePrinter
from pysisyphus.helpers import check_for_stop_sign
from pysisyphus.optimizers.closures import lbfgs_closure

np.set_printoptions(suppress=True, precision=4)


DimerCycle = namedtuple("DimerCycle",
                        "org_coords trial_coords rot_coords trans_coords f0 f0_mod",
)


def make_unit_vec(vec1, vec2):
    """Return unit vector pointing from vec2 to vec1."""
    diff = vec1 - vec2
    return diff / np.linalg.norm(diff)


def perpendicular_force(force, vec):
    return force - force.dot(vec)*vec


def rotate_R1(coords, rad, N, theta, dR):
    """Only valid for rotation of R1!"""
    return coords + (N*np.cos(rad) + theta*np.sin(rad)) * dR


def get_curvature(f1, f2, N, dR):
    return (f2 - f1).dot(N) / (2*dR)


def get_f_mod(f, N, C):
    f_mod = -f.dot(N)*N
    if C < 0:
        f_mod = f + 2*f_mod
    return f_mod


def get_geom_getter(ref_geom, calc_setter):
    def geom_from_coords(coords):
        new_geom = ref_geom.copy()
        new_geom.coords = coords
        new_geom.set_calculator(calc_setter())
        return new_geom
    return geom_from_coords


def make_theta(f1, f2, N):
    f1_perp = perpendicular_force(f1, N)
    f2_perp = perpendicular_force(f2, N)
    f_perp = f1_perp - f2_perp
    theta = f_perp / np.linalg.norm(f_perp)
    return theta


def plot_modes(trans_mode, N, N_trial, N_rot, coords0):
    fig, ax = plt.subplots()

    x, y = coords0[:2]
    trans_x, trans_y = trans_mode[:2]
    xt = np.array((x, x+trans_x))
    yt = np.array((y, y+trans_y))
    trans_line = mpl.lines.Line2D(xdata=xt, ydata=yt, label="Trans mode")
    lbls = "N N_trial N_rot".split()
    for i, (N_, l) in enumerate(zip((N, N_trial, N_rot), lbls)):
        xn, yn = N_[:2]
        xn_data = np.array((x, x+xn))
        yn_data = np.array((y, y+yn))
        n_line = mpl.lines.Line2D(xn_data, yn_data, label=l, color=f"C{i}")
        ax.add_line(n_line)

    ax.scatter(x, y)
    ax.add_line(trans_line)
    ax.legend()
    ax.set_xlim(-0.5, 2)
    ax.set_ylim(0, 3)

    plt.show()


def dimer_method(geoms, calc_getter, N_init=None,
                 max_step=0.1, max_cycles=50,
                 trial_angle=5, angle_thresh=0.5, dR_base=0.01,
                 ana_2dpot=False):
    """Dimer method using steepest descent for rotation and translation.

    See
        # Original paper
        [1] https://doi.org/10.1063/1.480097
        # Improved dimer
        [2] https://doi.org/10.1063/1.2104507
        # Several trial rotations
        [3] https://doi.org/10.1063/1.1809574
        # Superlinear dimer
        [4] https://doi.org/10.1063/1.2815812

        To add:
            Comparing curvatures and adding π/2 if appropriate.

        Default parameters from [1]
            max_step = 0.1 bohr
            dR_base = = 0.01 bohr
    """
    # Parameters
    angle_thresh_rad = np.deg2rad(angle_thresh)

    header = "Cycle Curvature max(f0) rms(f0)".split()
    col_fmts = "int float float float".split()
    table = TablePrinter(header, col_fmts)

    geom_getter = get_geom_getter(geoms[0], calc_getter)

    assert len(geoms) in (1, 2), "geoms argument must be of length 1 or 2!"
    if len(geoms) == 2:
        geom1, geom2 = geoms
        dR = np.linalg.norm(geom1.coords - geom2.coords) / 2
        # N points from geom2 to geom1
        N = make_unit_vec(geom1.coords, geom2.coords)
        coords0 = geom2.coords + dR*N
        geom0 = geom_getter(coords0)
    # This handles cases where only one geometry is supplied. We use the
    # geometry as midpoint, select a random dimer direction and derive
    # geom1 and geom2 from it.
    else:
        geom0 = geoms[0]
        # Assign random unit vector and use dR_base to for dR
        coord_size = geom0.coords.size
        N = N_init
        if N is None:
            N = np.random.rand(coord_size)
        if ana_2dpot:
            N[2] = 0
        N /= np.linalg.norm(N)
        coords1 = geom0.coords + dR_base*N
        geom1 = geom_getter(coords1)
        coords2 = geom0.coords - dR_base*N
        geom2 = geom_getter(coords2)
        dR = dR_base

    dimer_cycles = list()

    print("Using N:", N)
    def f0_mod_getter(coords, N, C):
        results = geom0.get_energy_and_forces_at(coords)
        forces = results["forces"]
        return get_f_mod(forces, N, C)

    def restrict_step(step, max_step=max_step):
        step_max = np.abs(step).max()
        if step_max > max_step:
            factor = max_step / step_max
            step *= factor
        return step

    table.print_header()
    for i in range(max_cycles):
        f0 = geom0.forces
        f1 = geom1.forces
        f2 = 2*f0 - f1

        coords0 = geom0.coords
        coords1 = geom1.coords
        coords2 = geom2.coords

        # Get rotated endpoint geometries. The rotation takes place in a plane
        # spanned by N and theta. Theta is a unit vector perpendicular to N that
        # can be formed from the perpendicular components of the forces at the
        # endpoints.
        theta = make_theta(f1, f2, N)

        C = get_curvature(f1, f2, N, dR)
        # Derivative of the curvature, Eq. (29) in [2]
        # (f2 - f1) or -(f1 - f2)
        dC = 2*(f0-f1).dot(theta)/dR
        trial_rad = -0.5*np.arctan2(dC, 2*abs(C))

        # Trial rotation for finite difference calculation of rotational force
        # and rotational curvature.
        coords1_star = rotate_R1(coords0, trial_rad, N, theta, dR)
        f1_star = geom1.get_energy_and_forces_at(coords1_star)["forces"]
        f2_star = 2*f0 - f1_star
        N_star = make_unit_vec(coords1_star, coords0)
        coords2_star = coords0 - N_star*dR

        trial_coords = np.array((coords1_star, coords0, coords2_star))
        C_star = get_curvature(f1_star, f2_star, N_star, dR)
        theta_star = make_theta(f1_star, f2_star, N_star)

        b1 = 0.5 * dC
        a1 = (C - C_star + b1*np.sin(2*trial_rad)) / (1-np.cos(2*trial_rad))
        a0 = 2 * (C - a1)

        rad_min = 0.5 * np.arctan(b1/a1)
        C_min = a0/2 + a1*np.cos(2*rad_min) + b1*np.sin(2*rad_min)
        if C_min > C:
            rad_min += np.deg2rad(90)

        if np.abs(rad_min) > angle_thresh_rad:
            coords1_rot = rotate_R1(coords0, rad_min, N, theta, dR)
            N_rot = make_unit_vec(coords1_rot, coords0)
            coords2_rot = coords0 - N_rot*dR
            geom1.coords = coords1_rot
            geom2.coords = coords2_rot
        # Don't do rotation for small angles
        else:
            N_rot = N
            coords1_rot = coords1
            coords2_rot = coords2
            table.print("Rotation angle too small. Skipping rotation.")

        # Translation
        f0_mod = get_f_mod(f0, N_rot, C_min)

        if i == 0:
            trans_lbfgs = lbfgs_closure(f0_mod, f0_mod_getter,
                                        restrict_step=restrict_step)
        N_trans = f0_mod.copy()
        N_trans /= np.linalg.norm(N_trans)

        if C > 0:
            step = max_step*N_trans
        else:
            lbfgs_res = trans_lbfgs(coords0, N_rot, C_min)
            step = lbfgs_res[1]

        coords0_trans = coords0 + step
        coords1_trans = coords0_trans + dR*N_rot
        coords2_trans = coords0_trans - dR*N_rot

        # Save cycle information
        org_coords = np.array((coords1, coords0, coords2))
        rot_coords = np.array((coords1_rot, coords0, coords2_rot))
        trans_coords = np.array((coords1_trans, coords0_trans, coords2_trans))
        dc = DimerCycle(org_coords, trial_coords, rot_coords, trans_coords, f0, f0_mod)
        dimer_cycles.append(dc)

        # Update dimer coordinates for next cycle
        geom0.coords = coords0_trans
        geom1.coords = coords1_trans
        geom2.coords = coords2_trans
        N = N_rot

        f0_rms = np.sqrt(np.power(f0, 2).mean())
        f0_max = np.abs(f0).max()
        row_args = [i, C, f0_max, f0_rms]
        table.print_row(row_args)
        if f0_rms <= 1e-3 and f0_max <= 1.5e-3:
            table.print("Converged!")
            break

        if check_for_stop_sign():
            break
    return dimer_cycles


if __name__ == "__main__":
    run()
