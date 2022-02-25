from darts.properties.geothermal_properties import *
from darts.properties.geothermal_physics import *
from darts.models.physics.iapws.iapws_property import *
from darts.models.physics.iapws.custom_rock_property import *
from model_base import BaseModel
from darts.engines import value_vector

class Model(BaseModel):

    def __init__(self, n_points=100):
        # call base class constructor
        super().__init__()

        self.hcap = np.array(self.reservoir.mesh.heat_capacity, copy = False)
        self.conduction = np.array(self.reservoir.mesh.rock_cond, copy = False)
        self.hcap.fill(2200)
        self.conduction.fill(181.44)

        # Create property containers:
        self.property_container = property_container(phase_name=['water', 'steam', 'temperature', 'energy'], component_name=['H2O'])

        # Define properties in property_container (IAPWS is the default property package for Geothermal in DARTS)
        # Users can define their custom properties in custom_properties.py; several property examples are defined there.
        self.rock = [value_vector([1, 0, 273.15])]
        self.property_container.temperature = iapws_temperature_evaluator()
        self.property_container.water_enthalpy = iapws_water_enthalpy_evaluator()
        self.property_container.steam_enthalpy = iapws_steam_enthalpy_evaluator()
        self.property_container.water_saturation = iapws_water_saturation_evaluator()
        self.property_container.steam_saturation = iapws_steam_saturation_evaluator()
        self.property_container.water_relperm = iapws_water_relperm_evaluator()
        self.property_container.steam_relperm = iapws_steam_relperm_evaluator()
        self.property_container.water_density = iapws_water_density_evaluator()
        self.property_container.steam_density = iapws_steam_density_evaluator()
        self.property_container.water_viscosity = iapws_water_viscosity_evaluator()
        self.property_container.steam_viscosity = iapws_steam_viscosity_evaluator()
        self.property_container.rock_compaction = custom_rock_compaction_evaluator(self.rock)
        self.property_container.rock_energy = custom_rock_energy_evaluator(self.rock)

        self.physics = Geothermal_Custom(timer=self.timer, n_points=128, min_p=1,
                                  max_p=151, min_e=1000, max_e=55000, property_container=self.property_container)

        self.params.first_ts = 1e-2
        self.params.mult_ts = 2
        self.params.max_ts = 15

        # Newton tolerance is relatively high because of L2-norm for residual and well segments
        self.params.tolerance_newton = 1e-3
        self.params.tolerance_linear = 1e-3
        self.params.max_i_newton = 20
        self.params.max_i_linear = 30
        self.runtime = 3650
        self.inj = value_vector([0.999])

    def set_initial_conditions(self):
        self.physics.set_uniform_initial_conditions(self.reservoir.mesh, uniform_pressure=100,
                                                    uniform_temperature=273.15+75)

    def set_boundary_conditions(self):
        for i, w in enumerate(self.reservoir.wells):
            if i > 14:
                # w.control = self.physics.new_rate_water_inj(120, 308.15)
                # w.control = self.physics.new_mass_rate_water_inj(2000, 1801.5)
                w.control = self.physics.new_bhp_water_inj(120, 308)
            else:
                # w.control = self.physics.new_rate_water_prod(120)
                w.control = self.physics.new_bhp_prod(80)

    def plot_temp_layer_map(self, map_data, k, name):
        import plotly
        import plotly.graph_objs as go

        nxny = self.reservoir.nx * self.reservoir.ny
        data = [go.Heatmap(
            z=map_data[nxny * (k - 1): nxny * k].reshape(self.reservoir.ny, self.reservoir.nx))]  # .transpose()
        layout = go.Layout(title='%s, layer %d' % (name, k),
                           yaxis=dict(scaleratio=1, scaleanchor='x', title='X, block'),
                           xaxis=dict(title='Y, block'))
        fig = go.Figure(data=data, layout=layout)
        plotly.offline.plot(fig, filename='%s_%d_map.html' % (name, k))

    def plot_pres_layer_map(self, map_data, k, name):
        import plotly
        import plotly.graph_objs as go

        nxny = self.reservoir.nx * self.reservoir.ny
        data = [go.Heatmap(
            z=map_data[nxny * (k - 1): nxny * k].reshape(self.reservoir.ny, self.reservoir.nx))]  # .transpose())]
        layout = go.Layout(title='%s, layer %d' % (name, k),
                           yaxis=dict(scaleratio=1, scaleanchor='x', title='X, block'),
                           xaxis=dict(title='Y, block'))
        fig = go.Figure(data=data, layout=layout)
        plotly.offline.plot(fig, filename='%s_%d_map.html' % (name, k))

    def enthalpy_to_temperature(self, data):
        from darts.models.physics.iapws.iapws_property_vec import _Backward1_T_Ph_vec
        data_len = int(len(data) / 2)
        T = np.zeros(data_len)
        T[:] = _Backward1_T_Ph_vec(data[::2]/10, data[1::2]/18.015)
        return T
        
    def export_vtk_geo(self, file_name='data'):
        from darts.models.physics.iapws.iapws_property_vec import _Backward1_T_Ph_vec
        nb = self.reservoir.nb
        T = np.zeros(nb)
        X = np.array(self.physics.engine.X)
        T[:] = _Backward1_T_Ph_vec(X[:nb*2:2]/10, X[1:nb*2:2]/18.015)
        self.export_vtk(file_name=file_name, local_cell_data={'temp': T})