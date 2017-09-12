#!/usr/bin/env python3
import logging
import os
import pathlib

import numpy as np
import pytest

from pysisyphus.calculators.ORCA import ORCA
from pysisyphus.Geometry import Geometry
from pysisyphus.irc.GonzalesSchlegel import GonzalesSchlegel

from qchelper.geometry import parse_xyz_file

@pytest.mark.orca_irc
def test_hcn_iso_gs_irc():
    this_dir = pathlib.Path(os.path.dirname(os.path.realpath(__file__)))
    ts_xyz_fn = this_dir / "hcn_ts_standard.xyz"
    atoms, coords = parse_xyz_file(ts_xyz_fn)
    geometry = Geometry(atoms, coords.flatten())
    geometry.set_calculator(ORCA())
    hessian = geometry.hessian

    gs_irc = GonzalesSchlegel(geometry, max_cycles=4, max_step=0.2, forward=True)
    gs_irc.run()
    gs_irc.write_trj(this_dir)


if __name__ == "__main__":
    test_hcn_iso_gs_irc()
