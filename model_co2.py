from properties.comp_physics import CO2Brine
from properties.comp_properties import *
from model_base import BaseModel
from darts.engines import value_vector
from properties.comp_operator_python import properties_evaluator

class Model(BaseModel):

    def __init__(self, n_points=100):
        # call base class constructor
        super().__init__()

        """Physical properties"""
        # Create property containers:
        self.property_container = property_container(phase_name=['gas', 'Aq'], component_name=['CO2', 'H2O'])
        self.components = self.property_container.component_name
        self.phases = self.property_container.phase_name

        """ properties correlations """
        self.property_container.flash_ev = Flash(self.components)
        self.property_container.density_ev = dict([('Aq', DensityBrine()), ('gas', DensityVap())])
        self.property_container.viscosity_ev = dict([('Aq', ViscosityBrine()), ('gas', ViscosityVap())])
        self.property_container.rel_perm_ev = dict([('Aq', PhaseRelPerm("Aq")), ('gas', PhaseRelPerm("gas"))])
        self.property_container.rel_well_perm_ev = dict([('Aq', WellPhaseRelPerm("Aq")), ('gas', WellPhaseRelPerm("gas"))])

        """ Activate physics """
        self.physics = CO2Brine(self.timer, n_points=501, min_p=100, max_p=500, min_z=1e-10, max_z=1-1e-10,
                                property_container=self.property_container)

        self.params.first_ts = 1e-2
        self.params.mult_ts = 2
        self.params.max_ts = 5

        # Newton tolerance is relatively high because of L2-norm for residual and well segments
        self.params.tolerance_newton = 1e-3
        self.params.tolerance_linear = 1e-3
        self.params.max_i_newton = 20
        self.params.max_i_linear = 30
        self.runtime = 900
        self.inj = value_vector([0.999])

    def set_initial_conditions(self):
        """ initialize conditions for all scenarios"""
        self.physics.set_uniform_initial_conditions(self.reservoir.mesh, 220, [1e-8])

    def set_boundary_conditions(self):
        for i, w in enumerate(self.reservoir.wells):
            if i == 18:
                w.control = self.physics.new_rate_gas_inj(320, self.inj)
                # w.control = self.physics.new_bhp_inj(230, self.inj_stream)
            else:
                w.control = self.physics.new_bhp_prod(210)

    def export_vtk_comp(self, file_name='data'):
        properties_vector = properties_evaluator(self.property_container)
        nb = self.reservoir.nb
        X = np.array(self.physics.engine.X)
        pres = X[0:2*nb:2]
        zg = X[1:2*nb:2]
        Sg = np.zeros(nb)
        xg = np.zeros(nb)

        for i in range (nb):
            state = value_vector([pres[i], zg[i]])
            [Sg[i], xg[i], ] = properties_vector.evaluate(state)


        self.export_vtk(file_name=file_name, local_cell_data={'GasSat': Sg, 'xCO2': xg})

