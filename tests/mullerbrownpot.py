#!/usr/bin/env python3

import numpy as np

from AnimPlot import AnimPlot
from calculators.MullerBrownPot import MullerBrownPot
from cos.NEB import NEB
from cos.SimpleZTS import SimpleZTS
from Geometry import Geometry
from optimizers.SteepestDescent import SteepestDescent
from optimizers.NaiveSteepestDescent import NaiveSteepestDescent

CYCLES = 20
IMAGES = 7


def get_geoms():
    min_a = np.array((-0.558, 1.442, 0)) # Minimum A
    min_b = np.array((0.6215, 0.02838, 0)) # Minimum B
    min_c = np.array((-0.05, 0.467, 0)) # Minimum C
    saddle_a = np.array((-0.822, 0.624, 0)) # Saddle point A
    coords = (min_b, min_c, saddle_a, min_a)
    atoms = ("H", )
    geoms = [Geometry(atoms, c) for c in coords]
    return geoms


def run_cos_opt(cos_class, reparametrize=False):
    geoms = get_geoms()
    cos = cos_class(geoms)
    cos.interpolate(IMAGES)
    for img in cos.images[1:-1]:
        img.set_calculator(MullerBrownPot())

    #sd = NaiveSteepestDescent(cos,
    sd = SteepestDescent(cos,
                         max_cycles=CYCLES,
                         max_force_thresh=0.05,
                         rms_force_thresh=0.01,
                         alpha=-0.005)
    if reparametrize:
        sd.run(reparam=cos.reparametrize)
    else:
        sd.run()
    xlim = (-1.75, 1.25)
    ylim = (-0.5, 2.25)
    levels=(-150, 5, 20)
    ap = AnimPlot(MullerBrownPot(), sd, xlim=xlim, ylim=ylim, levels=levels)
    ap.animate()


if __name__ == "__main__":
    #run_cos_opt(NEB)
    run_cos_opt(SimpleZTS, reparametrize=True)