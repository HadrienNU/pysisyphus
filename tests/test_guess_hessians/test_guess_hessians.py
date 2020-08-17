#!/usr/bin/env python3

import pytest

from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.helpers import geom_from_library, do_final_hessian
from pysisyphus.optimizers.guess_hessians import ts_hessian
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.tsoptimizers.RSPRFOptimizer import RSPRFOptimizer
from pysisyphus.testing import using

import matplotlib.pyplot as plt
import numpy as np


@using("pyscf")
@pytest.mark.parametrize(
    "hessian_init, ref_cycle", [
        # ("calc", 27),  # Very slow
        # Converges to wrong minimum
        # pytest.param("unit", 17, marks=[using("pyscf"), pytest.mark.xfail]),
        ("fischer", 17),
        ("lindh", 18),
        ("simple", 27),
        ("swart", 20),
        pytest.param("xtb", 22, marks=[using("pyscf"), using("xtb")]),
        pytest.param("xtb1", 20, marks=[using("pyscf"), using("xtb")]),
    ]
)
def test_guess_hessians(hessian_init, ref_cycle):
    geom = geom_from_library("birkholz/vitamin_c.xyz", coord_type="redund")
    geom.set_calculator(PySCF(basis="321g", pal=2))

    print("@\tguess_hessian:", hessian_init)
    opt_kwargs = {
        "hessian_init": hessian_init,
    }
    opt = RFOptimizer(geom, **opt_kwargs)
    opt.run()

    print("@\tcur_cycle:", opt.cur_cycle)

    assert opt.is_converged
    assert opt.cur_cycle == ref_cycle
    assert geom.energy == pytest.approx(-677.21287468)


@using("pyscf")
def test_ts_hessian():
    H = np.diag((1, 0.5, 0.25))
    tsh = ts_hessian(H, (0, 1))
    w, v = np.linalg.eigh(tsh)
    np.testing.assert_allclose(w, (-0.445194, 0.070194, 0.25), atol=1e-7)


@using("pyscf")
def test_ts_hessian_opt():
    geom = geom_from_library("baker_ts/01_hcn.xyz", coord_type="redund",)
    geom.set_calculator(PySCF(basis="321g"))

    opt_kwargs = {
        "hessian_init": "fischer",
        "dump": True,
        "rx_coords": ((2, 1, 0), ),
        "thresh": "gau_tight",
    }
    opt = RSPRFOptimizer(geom, **opt_kwargs)
    opt.run()
    do_final_hessian(geom)

    assert opt.is_converged
    assert opt.cur_cycle == 11
    assert geom.energy == pytest.approx(-92.2460426792319)