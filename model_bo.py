from properties.black_oil_physics import BlackOil
from properties.black_oil_properties import *
from properties.black_oil_operator_python import Properties
from model_base import BaseModel
from darts.engines import value_vector, sim_params

class Model(BaseModel):

    def __init__(self, n_points=100):
        # call base class constructor
        super().__init__()

        self.pvt = 'bo_physics.in'
        self.property_container = PropertyContainer(phase_name=['gas', 'oil', 'water'], component_name=['gas', 'oil', 'water'],
                                                    pvt=self.pvt)

        # Define property evaluators based on custom properties
        self.property_container.pb_ev = BubblePointPres(self.pvt)
        self.property_container.rs_ev = Rs(self.pvt)
        self.property_container.xgo_ev = Xgo(self.pvt)
        self.property_container.density_ev = dict([('gas', DensityGas(self.pvt)), ('oil', DensityOil(self.pvt)),
                                                   ('water', DensityWater(self.pvt))])
        self.property_container.viscosity_ev = dict([('gas', VsicoGas(self.pvt)), ('oil', ViscoOil(self.pvt)),
                                                     ('water', ViscoWat(self.pvt))])
        self.property_container.sat_ev = dict([('gas', Gassat(self.pvt)), ('oil', Oilsat(self.pvt)),
                                                     ('water', Watsat(self.pvt))])
        # stone I model is used, then OW is water-oil system, and GO is gas-oil system
        self.property_container.rel_perm_ev = dict([('gas', GasRelPerm(self.pvt)), ('oil', OilRelPerm(self.pvt)),
                                                    ('water', WatRelPerm(self.pvt)), ('OW', Krow(self.pvt)),
                                                    ('GO', Krog(self.pvt))])
        self.property_container.rock_compress_ev = RockCompactionEvaluator(self.pvt)

        # create physics
        self.grav = 0
        self.physics = BlackOil(self.timer, n_points=500, min_p=0, max_p=450, min_z=1e-13,
                                property_container=self.property_container, grav=self.grav)

        self.params.first_ts = 1e-2
        self.params.mult_ts = 2
        self.params.max_ts = 15

        # Newton tolerance is relatively high because of L2-norm for residual and well segments
        self.params.tolerance_newton = 1e-2
        self.params.tolerance_linear = 1e-3
        self.params.max_i_newton = 20
        self.params.max_i_linear = 30
        self.params.newton_type = sim_params.newton_local_chop
        self.params.nonlinear_norm_type = sim_params.L1
        self.inj = value_vector([1e-8, 1e-8])

    def set_initial_conditions(self):
        self.physics.set_uniform_initial_conditions(self.reservoir.mesh, uniform_pressure=320,
                                                    uniform_composition=[0.007296, 0.4832])

    def set_boundary_conditions(self):
        for i, w in enumerate(self.reservoir.wells):
            if i > 14:
                w.control = self.physics.new_bhp_inj(400, self.inj)
            else:
                w.control = self.physics.new_rate_oil_prod(3000)
                w.constraint = self.physics.new_bhp_prod(70)

    def export_vtk_bo(self, file_name='data'):
        properties_vector = Properties(self.property_container)
        nb = self.reservoir.nb
        X = np.array(self.physics.engine.X)
        pres = X[0:3*nb:3]
        zg = X[1:3*nb:3]
        zo = X[2:3*nb:3]
        So = np.zeros(nb)
        Sg = np.zeros(nb)
        Sw = np.zeros(nb)
        pbub = np.zeros(nb)

        for i in range (nb):
            state = value_vector([pres[i], zg[i], zo[i]])
            [So[i], Sw[i], Sg[i], pbub[i]] = properties_vector.evaluate(state)


        self.export_vtk(file_name=file_name, local_cell_data={'OilSat': So, 'GasSat': Sg, 'WatSat': Sw})

