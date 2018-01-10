#!/usr/bin/env python3

import os
import re

import numpy as np

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.config import Config
from pysisyphus.constants import BOHR2ANG

from pysisyphus.xyzloader import make_xyz_str


class OpenMolcas(Calculator):

    def __init__(self, basis, inporb, roots, mdrlxroot,
                 supsym=None, **kwargs):
        super(OpenMolcas, self).__init__(**kwargs)

        self.basis = basis
        self.inporb = inporb
        self.roots = roots
        self.mdrlxroot = mdrlxroot
        self.supsym = self.build_supsym_str(supsym)
        self.cur_jobiph = ""
        self.prev_jobiph = ""

        self.inp_fn = "openmolcas.in"
        self.out_fn = "openmolcas.out"
        self.float_regex = "([\d\.\-]+)"

        self.openmolcas_input = """
        >> copy {inporb}  $Project.RasOrb
        &gateway
         coord
          {xyz_str}
         basis
          {basis}
         group
          nosym
        ricd
 
        &seward
         doanalytical

        &rasscf
         charge
          {charge}
         spin
          {mult}
         fileorb
          $Project.RasOrb
         thrs
          1.0e-6,1.0e-2,1.0e-2
         ciroot
          {roots} {roots} 1
         mdrlxroot
          {mdrlxroot}
         {supsym}

        >> copy $Project.JobIph $CurrDir/$Project.JobIph

        {rassi}

        &alaska
         pnew
        """

        self.parser_funcs = {
            "grad": self.parse_gradient,
        }

        self.base_cmd = Config["openmolcas"]["cmd"]

    def build_supsym_str(self, supsym):
        """Can handle only one subgroup for now."""
        if not supsym:
            return ""

        num_orbitals = len(supsym.split())
        return f"supsym\n1\n{num_orbitals} {supsym};"

    def build_rassi_str(self):
        if self.counter == 0:
            return ""
        else:
            return f"""
            >> copy $Project.JobIph JOB001
            >> copy {self.prev_jobiph} JOB002
            &rassi
             track
            """

    def prepare_coords(self, atoms, coords):
        coords = coords * BOHR2ANG
        return make_xyz_str(atoms, coords.reshape((-1, 3)))

    def prepare_input(self, atoms, coords):
        xyz_str = self.prepare_coords(atoms, coords)
        inp = self.openmolcas_input.format(
                                        inporb=self.inporb,
                                        xyz_str=xyz_str,
                                        basis=self.basis,
                                        charge=self.charge,
                                        mult=self.mult,
                                        roots=self.roots,
                                        mdrlxroot=self.mdrlxroot,
                                        supsym=self.supsym,
                                        rassi=self.build_rassi_str(),
        )
        return inp

    def get_forces(self, atoms, coords):
        self.logger.debug(f"using inporb: {self.inporb}")
        inp = self.prepare_input(atoms, coords)
        add_args = ("-clean", "-oe", self.out_fn)
        env = os.environ.copy()
        env["MOLCAS_PROJECT"] = f"{self.name}_{self.counter}"
        results = self.run(inp, calc="grad", add_args=add_args, env=env)
        return results

    def keep(self, path):
        kept_fns = super().keep(path, ("RasOrb", "out", "in", "JobIph"))
        self.inporb = kept_fns["RasOrb"]
        # Keep references to the current and the last .JobIph file
        # to be used in &rassi to track our root in a state average
        # calculation.
        # In the first iteration self.cur_jobiph isn't set yet
        if self.counter == 0:
            self.prev_jobiph = kept_fns["JobIph"]
        else:
            self.prev_jobiph = self.cur_jobiph
        self.cur_jobiph = kept_fns["JobIph"]
        self.log(f"current JobIph {self.cur_jobiph}")

    def parse_energies(self, text):
        # Energy of root for which gradient was computed
        energy_regex = "RASSCF state energy =\s*" + self.float_regex
        energy = float(re.search(energy_regex, text).groups()[0])

        # All state average energies
        root_re = "RASSCF root number.+Total energy.+?" + self.float_regex
        matches = re.findall(root_re, text)
        sa_energies = np.array(matches, dtype=np.float)

        return energy, sa_energies

    def parse_gradient(self, path):
        results = {}
        gradient_fn = os.path.join(path, self.out_fn)
        with open(gradient_fn) as handle:
            text = handle.read()

        # Search for the block containing the gradient table
        regex = "Molecular gradients(.+?)--- Stop Module:\s*alaska"
        floats = [self.float_regex for i in range(3)]
        line_regex = "([A-Z\d]+)\s*" + "\s*".join(floats)

        mobj = re.search(regex, text, re.DOTALL)
        gradient = list()
        for line in mobj.groups()[0].split("\n"):
            # Now look for the lines containing the gradient
            mobj = re.match(line_regex, line.strip())
            if not mobj:
                continue
            # Discard first column (atom+number)
            gradient.append(mobj.groups()[1:])
        gradient = np.array(gradient, dtype=np.float).flatten()

        if self.counter > 0:
            self.parse_rassi_track(path)

        energy, sa_energies = self.parse_energies(text)

        results["energy"] = energy
        results["sa_energies"] = sa_energies
        results["forces"] = -gradient

        return results

    def parse_rassi_track(self, path):
        gradient_fn = path / self.out_fn
        with open(gradient_fn) as handle:
            text = handle.read()
        track_re = "Initial root:\s*(\d+)\s*Overlaps with current " \
                   "states:(.+)New root:\s*(\d+)"
        #overlap_re = "OVERLAP MATRIX FOR THE ORIGINAL STATES:(.+?)##"
        mobj = re.search(track_re, text, re.DOTALL)

        initial_root, overlaps, new_root = mobj.groups()
        overlaps = np.array(overlaps.strip().split(), dtype=float).reshape(-1, 2)
        # Filters for overlaps > 10% (0.1**2 ~ 0.31622)
        thresh = 0.1
        inds = np.where(np.abs(overlaps[:,1]) > thresh**0.5)
        ov_perc_str = ", ".join([f"{nr:.0f}: {ov**2:.2%}"
                                 for nr, ov in overlaps[inds]])
        self.log(f"Overlaps between previous root {initial_root} and "
                 f"new roots bigger {thresh:.0%}:  {ov_perc_str}. Will "
                 f"use root {new_root} for the following gradient calculation.")
        if new_root != initial_root:
            self.log("Found a root flip!")
        self.mdrlxroot = new_root

    def __str__(self):
        return "OpenMolcas calculator"


if __name__ == "__main__":
    from pysisyphus.helpers import geom_from_library
    fileorb = "/scratch/test/ommin/excrp.es_opt.RasOrb"
    basis = "6-31G*"
    roots = 5
    rlxroot = 5
    om = OpenMolcas(basis, fileorb, roots, rlxroot)
    geom = geom_from_library("dieniminium_cation_s1_opt.xyz")
    geom.set_calculator(om)
    #print(geom.forces)
    #om.parse_gradient("/scratch/test/satest")
    from pathlib import Path
    p = Path("/scratch/track_test/parse")
    om.parse_rassi_track(p)