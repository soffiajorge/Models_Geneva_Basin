from darts.models.reservoirs.struct_reservoir import StructReservoir
from darts.models.physics.dead_oil import DeadOil
from model_base import BaseModel
from darts.engines import value_vector, sim_params
import numpy as np
from darts.tools.keyword_file_tools import load_single_keyword
import os

class Model(BaseModel):

    def __init__(self, n_points=100):
        # call base class constructor
        super().__init__()

        self.physics = DeadOil(self.timer, 'physics.in', n_points, 0, 500, 1e-8)    #1e4
        self.params.first_ts = 1e-4
        self.params.mult_ts = 6
        self.params.max_ts = 15

        # Newton tolerance is relatively high because of L2-norm for residual and well segments
        self.params.tolerance_newton = 1e-3
        self.params.tolerance_linear = 1e-5
        self.params.max_i_newton = 20
        self.params.max_i_linear = 30
        self.runtime = 900

    def set_initial_conditions(self):
        self.physics.set_uniform_initial_conditions(self.reservoir.mesh, uniform_pressure=250,
                                                    uniform_composition=[0.2357])

    def set_boundary_conditions(self):
        for i, w in enumerate(self.reservoir.wells):
            if i >= 14:
                w.control = self.physics.new_rate_water_inj(2000)
                w.constraint = self.physics.new_bhp_water_inj(300)
            else:
                w.control = self.physics.new_bhp_prod(150)

