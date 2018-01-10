#!/usr/bin/env python3

import numpy as np

from pysisyphus.helpers import fit_rigid
from pysisyphus.optimizers.Optimizer import Optimizer


class QuickMin(Optimizer):

    def __init__(self, geometry, dt=0.35, **kwargs):
        super(QuickMin, self).__init__(geometry, **kwargs)

        self.dt = dt

    def prepare_opt(self):
        self.velocities = [np.zeros_like(self.geometry.coords), ]

    def optimize(self):
        if self.align and self.is_cos:
            (self.velocities[-1], ) = fit_rigid(self.geometry,
                                                (self.velocities[-1], ))

        prev_velocities = self.velocities[-1]
        cur_forces = self.geometry.forces
        self.forces.append(cur_forces)
        self.energies.append(self.geometry.energy)

        if self.cur_cycle == 0:
            tmp_velocities = np.zeros_like(prev_velocities)
        else:
            overlap = prev_velocities.dot(cur_forces)
            if overlap > 0:
                tmp_velocities = (overlap * cur_forces
                                  / cur_forces.dot(cur_forces))
            else:
                tmp_velocities = np.zeros_like(prev_velocities)
                self.log("resetted velocities")

        accelerations = cur_forces / self.geometry.masses_rep
        cur_velocities = tmp_velocities + self.dt*accelerations
        steps = cur_velocities*self.dt + 1/2*accelerations*self.dt**2
        steps = self.scale_by_max_step(steps)
        self.velocities.append(cur_velocities)
        velo_norm = np.linalg.norm(cur_velocities)
        acc_norm = np.linalg.norm(accelerations)
        self.log(f"norm(v) = {velo_norm:.4f}, norm(a) = {acc_norm:.4f}")

        return steps