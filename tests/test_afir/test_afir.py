#!/usr/bin/env python3

from pysisyphus.calculators.AFIR import AFIR
from pysisyphus.calculators.XTB import XTB
from pysisyphus.helpers import geom_from_library
from pysisyphus.optimizers.RFOptimizer import RFOptimizer


def test_afir():
    # geom = geom_from_library("ohch3f_anion_cs.xyz", coord_type="redund")
    geom = geom_from_library("ohch3f_anion_cs.xyz")
    fragment_indices = ([0, 1, 2, 3, 4], [5, 6])
    calc = XTB(charge=-1)
    gamma = 100 * 3.8088e-4
    afir = AFIR(calc, geom.atoms, fragment_indices, gamma)
    geom.set_calculator(afir)

    opt = RFOptimizer(geom, dump=True, trust_max=.3)
    opt.run()
    assert opt.is_converged


if __name__ == "__main__":
    test_afir()
