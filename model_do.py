from properties.dead_oil_physics import DeadOil
from properties.dead_oil_properties import *
from model_base import BaseModel
from darts.engines import value_vector

class Model(BaseModel):
    def __init__(self, n_points=100):
        # call base class constructor
        super().__init__()

        property_container = PropertyContainer(phase_name=['water', 'oil'], component_name=['w', 'o'])

        # Define property evaluators based on custom properties
        property_container.density_ev = dict([('water', DensityWater()), ('oil', DensityOil())])
        property_container.viscosity_ev = dict([('water', ViscosityWater()), ('oil', ViscosityOil())])
        property_container.watersat_ev = WaterSaturation()
        property_container.rel_perm_ev = dict([('water', PhaseRelPerm("water")), ('oil', PhaseRelPerm("oil"))])
        property_container.capillary_ev = CapillaryPressure()
        property_container.rock_compress_ev = RockCompactionEvaluator()

        # create physics
        self.grav = 0
        self.physics = DeadOil(self.timer, n_points=400, min_p=0, max_p=1000, min_z=1e-13,
                               property_container=property_container, grav=self.grav)
        self.params.first_ts = 1e-2
        self.params.mult_ts = 2
        self.params.max_ts = 15

        # Newton tolerance is relatively high because of L2-norm for residual and well segments
        self.params.tolerance_newton = 1e-3
        self.params.tolerance_linear = 1e-3
        self.params.max_i_newton = 20
        self.params.max_i_linear = 30
        self.runtime = 900
        self.inj = value_vector([0.999])

    def set_initial_conditions(self):
        self.physics.set_uniform_initial_conditions(self.reservoir.mesh, uniform_pressure=250,
                                                    uniform_composition=[0.2357])

    def set_boundary_conditions(self):
        for i, w in enumerate(self.reservoir.wells):
            if i >= 14:
                w.control = self.physics.new_rate_water_inj(2000, self.inj)
                w.constraint = self.physics.new_bhp_water_inj(300, self.inj)
            else:
                w.control = self.physics.new_bhp_prod(150)

